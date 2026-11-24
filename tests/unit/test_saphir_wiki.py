"""saphir_wiki.py: from the text of a wiki day page to roof times, injections and checks.  Everything here works
on text written in the test; nothing is fetched."""
from types import SimpleNamespace as NS

import pytest

from juice.core import injection_calc as ic
from juice.importers import saphir_wiki as w

T, P = 300.0, 1000.0
CH = NS(volume_m3=ic.R_GAS * T / 1000.0, syringe_T_K=300.0, mfc_std_T_K=300.0, mfc_std_P_hPa=1000.0)   # 100 mol of air


# ── quantities in a line ─────────────────────────────────────────────────────

@pytest.mark.parametrize("text, minutes", [
    ("50 sec", 50 / 60), ("50 seconds at 5 slpm", 50 / 60), ("5 min", 5.0), ("3 minutes", 3.0), ("1,5 min", 1.5),
    ("OH reactivity 25 s-1", None), ("10 s -1", None), ("20 µl isoprene", None), ("3 sccm", None),
])
def test_duration(text, minutes):
    assert w.parse_duration_min(text) == (pytest.approx(minutes) if minutes is not None else None)


def amounts(text):
    return {k: v for k, v in w.parse_amounts(text).items() if v is not None}


@pytest.mark.parametrize("text, found", [
    ("20 µl isoprene by syringe", {"ul": 20.0}),
    ("500 nl + 1 µl", {"ul": 1.5}),
    ("2x3.2 ml propene (1%)", {"ml": 6.4, "bottle_ppm": 10000.0}),
    ("7,5 ppb NO", {"ppb": 7.5}),
    ("0.5 ppm O3", {"ppb": 500.0}),                                       # a stated mixing ratio in ppm
    ("100 ppm bottle, 3 ml", {"ml": 3.0, "bottle_ppm": 100.0}),            # ppm of a bottle, not of the chamber
    ("NO (3835ppm) 50 s at 5 slpm", {"sccm": 5000.0, "duration_s": 50.0, "bottle_ppm": 3835.0}),
    ("FC 0.5l 62.5% for 2 min", {"sccm": 312.5, "duration_s": 120.0}),     # flow controller: full scale × setting
    ("Phenolic Mix 43ul for 234ppbc 25 s-1", {"ul": 43.0, "ppbc": 234.0}),
    ("roof open", {}),
])
def test_amounts_in_a_line(text, found):
    assert amounts(text) == pytest.approx(found)


# ── the action of a line without a code ──────────────────────────────────────

@pytest.mark.parametrize("text, code", [
    ("H.F. CHANGE: Roof open", "R"), ("R.W. CHANGE: roof closed", "r"), ("Dach auf", "R"), ("Dach zu", "r"),
    ("Ö.E. ZUSTAND: roof open", ""),                                       # a state report, not an action
    ("humidification started", "H"), ("Humidification off", "h"), ("seed aerosol on", "S"), ("seed stop", "s"),
    ("CHANGE: Injection of 20 µl isoprene", "X"), ("added 40 ppb ozone", "O"), ("NO addition 2 ppb", "M"),
    ("NO2 5 ppb", "D"), ("CO 500 ppb added", "c"), ("CO2 injection", "C"),
    ("flushing started", ""), ("fan on", ""), ("roof", ""),
])
def test_action_is_read_from_the_words(text, code):
    assert w.infer_code(text) == code


def test_the_word_no_is_not_nitric_oxide():
    """Bug of 24 Nov 2026: 'no injection today' was read as an NO injection and imported with the default 5 ppb."""
    assert w.infer_code("no injection today") == ""
    assert w.infer_code("CHANGE: No addition of VOC, instrument problem") == ""
    assert w.infer_code("injection of NO, 2 ppb") == "M" and w.infer_code("nitric oxide added") == "M"


def test_action_names_for_the_import_dialog():
    names = {code: w.WikiEvent("10:00", code, "").action for code in "RrHhXOMcSs$"}
    assert names["R"] == "roof open" and names["r"] == "roof close" and names["$"] == "note"
    assert names["O"] == "Ozone injection" and names["X"] == "VOC injection"
    assert names["h"].endswith("— stop") and names["H"] == names["h"][:-len(" — stop")]


# ── the action log of a page ─────────────────────────────────────────────────

PAGE = """<html><head><style>p {color: red}</style></head><body><h1>2024-07-09</h1>
<p>08:05 R H.F. CHANGE: Roof&nbsp;open</p>
<p>9:30:15 X  R.W.   CHANGE: 20 µl isoprene</p>
<p>10:00 H.F. CHANGE: added 40 ppb ozone</p>
<p>10:30 Ö.E. ZUSTAND: everything fine</p>
<p>11:00 25.3 21.0 1013</p>
<p>25:00 this is not a time of day</p>
<p>Plan for tomorrow: 09:00 start</p>
<p>14:30 r H.F. CHANGE: Roof closed</p></body></html>"""


def test_coded_lines_keep_their_code_and_uncoded_ones_get_one():
    events = w.parse_action_log(PAGE)
    assert [(e.time, e.code) for e in events] == [("08:05", "R"), ("09:30", "X"), ("10:00", "O"), ("10:30", "-"), ("14:30", "r")]
    assert events[0].description == "H.F. CHANGE: Roof open"
    assert events[1].description == "R.W. CHANGE: 20 µl isoprene"         # several blanks become one
    assert events[2].line == "10:00 H.F. CHANGE: added 40 ppb ozone"       # the line as it was read


def test_raw_wiki_text_is_read_the_same_way():
    raw = "08:05 R H.F. CHANGE: Roof open\n14:30 r H.F. CHANGE: Roof closed\n"
    assert [(e.time, e.code) for e in w.parse_action_log(raw)] == [("08:05", "R"), ("14:30", "r")]
    assert w.strip_html("no tags & nothing to do") == "no tags & nothing to do"


# ── the amount of one compound ───────────────────────────────────────────────

def compound(mcm="X", mw=50.0, density=0.5):
    return ic.ResolvedCompound("x", mcm, mw, density, "mcm")


def test_liquid_volume_becomes_ppb():
    ppb, delivery, text, alerts = w.amount_for_compound("10 µl by syringe", compound(), T, P, CH)
    assert ppb == pytest.approx(1000.0) and alerts == []                   # 1e-4 mol in 100 mol of air
    assert delivery == {"mode": "liquid", "vol_ul": "10", "mw": "50", "density": "0.5"}


def test_stated_ppb_wins_and_a_calculation_that_disagrees_is_flagged():
    ppb, delivery, text, alerts = w.amount_for_compound("10 µl for 1000 ppb", compound(), T, P, CH)
    assert ppb == 1000.0 and delivery["mode"] == "ppb" and alerts == [] and "as stated in the wiki" in text
    ppb, _delivery, _text, alerts = w.amount_for_compound("10 µl for 500 ppb", compound(), T, P, CH)
    assert ppb == 500.0 and "differ by more than 15 %" in alerts[0]


def test_gas_volume_without_bottle_concentration_is_taken_as_pure_gas():
    ppb, delivery, _text, alerts = w.amount_for_compound("3 ml", compound(density=None), T, P, CH)
    assert ppb == pytest.approx(3e-6 / CH.volume_m3 * 1e9) and delivery["bottle_unit"] == "pure"
    assert alerts == ["gas volume given without bottle concentration: pure gas assumed"]
    ppb, delivery, _text, alerts = w.amount_for_compound("3 ml (1%)", compound(density=None), T, P, CH)
    assert ppb == pytest.approx(3e-6 / CH.volume_m3 * 1e9 / 100.0) and alerts == []


def test_flow_and_time_need_the_bottle_concentration():
    ppb, delivery, _text, alerts = w.amount_for_compound("60 sccm for 60 s", compound(), T, P, CH)
    assert ppb is None and delivery["mode"] == "mfc" and "no bottle concentration" in alerts[0]
    ppb, _delivery, _text, alerts = w.amount_for_compound("NO (1000ppm) 60 sccm for 60 s", compound(), T, P, CH)
    assert ppb == pytest.approx(60e-6 / CH.volume_m3 * 1e9 * 1e-3)


def test_volume_of_a_compound_without_molar_mass_and_a_line_without_amount():
    ppb, _delivery, _text, alerts = w.amount_for_compound("10 µl", compound(mw=None, density=None), T, P, CH)
    assert ppb is None and "molar mass/density unknown" in alerts[0] and alerts[-1] == "no amount found"
    ppb, _delivery, text, alerts = w.amount_for_compound("some of it", compound(), T, P, CH)
    assert ppb is None and "enter ppb" in text and alerts == ["no amount found"]


def test_names_for_mixtures_without_a_recipe():
    assert w.mixture_label("R.W. CHANGE: Inject Phenolic + Furanoid Mix 43ul for 234ppbc") == "Phenolic_Furanoid_Mix"
    assert w.mixture_label("CHANGE: 43 ul") == ""
    assert w.identify_mixture_word("diesel exhaust added") == "DIESEL" and w.identify_mixture_word("VCP solution") == "VCP"
    assert w.identify_mixture_word("mixing fan on") == ""                  # 'mix' only as a word of its own


# ── from the log to the experiment ───────────────────────────────────────────

def imported(lines):
    events = w.parse_action_log("\n".join(lines))
    return w.events_to_import(events, "2024-07-09", T, P, "SAPHIR", mixtures={})


def checks(imp, part):
    return [text for _level, text in imp.checks if part in text]


def test_roof_times_are_paired_in_order():
    imp = imported(["08:00 R open", "10:00 r closed", "11:00 R open", "13:00 r closed", "09:00 M 2 ppb NO"])
    assert imp.roof_open == [("08:00", "10:00"), ("11:00", "13:00")] and not checks(imp, "roof")


def test_roof_entries_without_their_partner():
    imp = imported(["07:00 r closed", "08:00 R open", "09:00 M 2 ppb NO"])
    assert imp.roof_open == [("08:00", "23:59")]
    assert checks(imp, "07:00: roof closed without a preceding open") and checks(imp, "08:00: roof opened without a close")


def test_humidification_from_start_to_stop_with_its_impurities():
    imp = imported(["06:00 H start", "06:40 h stop"])
    (h2o,) = imp.injections
    assert (h2o.name, h2o.action, h2o.start, h2o.duration_min, h2o.column) == ("H2O", "humidification", "06:00", 40.0, "QH2O")
    assert [p for p, _amount in h2o.extra_channels] == ["CO", "NO2+XNOy", "NO+XNOy"]


def test_stop_logged_with_the_start_code_and_starts_without_stop():
    imp = imported(["06:00 H humidification", "06:30 H humidification stop"])
    assert [(e.start, e.duration_min) for e in imp.injections] == [("06:00", 30.0)]
    imp = imported(["06:00 H first", "07:00 H second"])
    assert [(e.name, e.start, e.duration_min) for e in imp.injections] == [("H2O", "06:00", 45.0), ("H2O_2", "07:00", 45.0)]
    assert checks(imp, "started twice") and checks(imp, "07:00: Humidification") and imp.injections[1].column == "QH2O_2"
    assert checks(imported(["06:00 h stop"]), "stopped without a start")


def test_inorganic_injections_take_the_stated_amount_or_the_default():
    imp = imported(["09:00 M 2 ppb NO", "09:30 O ozone added", "10:00 M 3 ppb NO for 5 min"])
    no, o3, no_2 = imp.injections
    assert (no.name, no.products, no.total_ppb, no.column, no.duration_min) == ("NO", "NO", "2", "QNO", 1.0)
    assert (o3.name, o3.total_ppb) == ("O3", "50") and checks(imp, "09:30: Ozone injection") and checks(imp, "default 50 ppb used")
    assert (no_2.name, no_2.column, no_2.switch, no_2.total_ppb, no_2.duration_min) == ("NO_2", "QNO_2", "qNO_2", "3", 5.0)


def test_voc_injection_becomes_its_mcm_species_and_a_precursor():
    imp = imported(["10:00 X 20 µl isoprene by syringe"])
    (voc,) = imp.injections
    assert (voc.name, voc.products, voc.action, voc.column) == ("C5H8", "C5H8", "voc", "Q_C5H8")
    assert imp.precursors == ["C5H8"] and voc.delivery["mode"] == "liquid" and float(voc.total_ppb) > 0
    assert voc.comment.isascii()                                           # the EASY file takes ASCII only


def test_compound_outside_the_mcm_and_text_that_nobody_knows_are_recorded_only():
    imp = imported(["10:00 X 5 µl texanol", "10:30 X the usual stuff"])
    texanol, unknown = imp.injections
    assert (texanol.name, texanol.products, texanol.action) == ("texanol", "", "voc") and checks(imp, "'texanol' has no MCM species")
    assert (unknown.name, unknown.products, unknown.action) == ("VOC", "", "mixture") and checks(imp, "not recognised")
    assert imp.precursors == []


def test_mixture_stated_in_ppbc_without_a_recipe():
    imp = imported(["10:00 X R.W. CHANGE: Inject Phenolic Mix 43ul for 234ppbc 25 s-1"])
    (mix,) = imp.injections
    assert (mix.name, mix.action, mix.products) == ("Phenolic_Mix", "mixture", "")
    assert "234 ppbC stated; OH reactivity 25 s-1" in mix.comment and checks(imp, "has no recipe file")


def test_mixture_with_a_recipe_is_expanded():
    recipe = ic.Mixture("VCP2", aliases=["VCP"], components=[ic.MixtureComponent("limonene", 10.0, "ul", "LIMONENE", "", 100.0, 1.0),
                                                              ic.MixtureComponent("ethanol", 30.0, "ul", "C2H5OH", "", 100.0, 1.0)])
    events = w.parse_action_log("10:00 X VCP solution (20 ul)")
    imp = w.events_to_import(events, "2024-07-09", T, P, "SAPHIR", mixtures={"VCP2": recipe})
    (mix,) = imp.injections
    assert mix.action == "mixture" and mix.products == "0.75*C2H5OH+0.25*LIMONENE" and mix.name == "VCP2"
    assert imp.precursors == ["LIMONENE", "C2H5OH"] and mix.delivery["injected_ul"] == "20"


def test_notes_unticked_lines_and_the_order_of_the_result():
    events = w.parse_action_log("12:00 M 2 ppb NO\n09:00 $ instrument check\n10:00 M 1 ppb NO\n11:00 X 5 µl isoprene")
    events[3].use = False                                                  # unticked in the import dialog
    imp = w.events_to_import(events, mixtures={})
    assert [e.start for e in imp.injections] == ["10:00", "12:00"] and imp.notes == ["09:00 $ instrument check"]
    assert [e.name for e in imp.injections] == ["NO_2", "NO"]              # named in page order, listed by time


def test_a_page_without_injections_says_so():
    assert checks(imported(["08:00 R open", "12:00 r closed"]), "no injection entries found")


# ── fetching ─────────────────────────────────────────────────────────────────

class FakeSession:
    """Answers like the wiki: 403 until logged in with the right password."""

    def __init__(self, page_status=200, password="secret"):
        self.page_status, self.password, self.logged_in, self.posts = page_status, password, False, []

    def get(self, url, params=None, timeout=None, verify=None):
        if self.page_status == 200 and not self.logged_in:
            return NS(status_code=403, content=b"You are not allowed to view this page.", text="")
        body = "08:00 R open" if self.page_status == 200 else ""
        return NS(status_code=self.page_status, content=body.encode(), text=body)

    def post(self, url, params=None, data=None, timeout=None, verify=None):
        self.posts.append(data)
        self.logged_in = data["password"] == self.password


def fetch(monkeypatch, session, user="", password=""):
    requests = pytest.importorskip("requests")
    monkeypatch.setattr(requests, "Session", lambda: session)
    monkeypatch.delenv(w.ENV_USER, raising=False)
    monkeypatch.delenv(w.ENV_PASS, raising=False)
    return w.fetch_wiki_text("2024-07-09", user, password)


def test_page_is_fetched_after_a_login(monkeypatch):
    session = FakeSession()
    assert fetch(monkeypatch, session, "me", "secret") == "08:00 R open"
    assert session.posts[0]["name"] == "me"


def test_page_that_needs_a_login_says_what_to_do(monkeypatch):
    with pytest.raises(PermissionError, match="Enter your wiki user name and password"):
        fetch(monkeypatch, FakeSession())
    with pytest.raises(PermissionError, match="login with user 'me' was not accepted"):
        fetch(monkeypatch, FakeSession(), "me", "wrong")


def test_missing_page_and_server_error(monkeypatch):
    with pytest.raises(FileNotFoundError, match="No wiki page for 2024-07-09"):
        fetch(monkeypatch, FakeSession(page_status=404))
    with pytest.raises(RuntimeError, match="HTTP 500"):
        fetch(monkeypatch, FakeSession(page_status=500))


def test_login_falls_back_to_the_environment(monkeypatch):
    monkeypatch.setenv(w.ENV_USER, " envuser ")
    monkeypatch.setenv(w.ENV_PASS, "envpass")
    assert w.wiki_credentials() == ("envuser", "envpass") and w.wiki_credentials("me", "pw") == ("me", "pw")
    assert w.wiki_url("2024-07-09").endswith("/saphirwiki/2024-07-09")
