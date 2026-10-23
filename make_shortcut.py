r"""
make_shortcut.py — put a JUICE icon on the Desktop (and optionally the Start menu).

    python make_shortcut.py              # Desktop\JUICE.lnk → JUICE.bat with the JUICE icon
    python make_shortcut.py --here       # JUICE.lnk inside the project folder instead
    python make_shortcut.py --startmenu  # also Start Menu\Programs\JUICE.lnk

Each click opens a new, independent JUICE window (JUICE.bat uses `start` + pythonw), so
several experiments can be worked on in parallel without a terminal.  Windows only.

A .lnk always stores absolute paths of the machine it was made on, so it is NOT synced
through git (JUICE.lnk is ignored); JUICE.bat is, and works anywhere.  After a clone or on
another PC run this script once to get the icon.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BAT = os.path.join(HERE, "JUICE.bat")
ICO = os.path.join(HERE, "JUICE.ico")


def make(link_path: str) -> None:
    ps = f'''
$s = (New-Object -ComObject WScript.Shell).CreateShortcut("{link_path}")
$s.TargetPath = "{BAT}"
$s.WorkingDirectory = "{HERE}"
$s.IconLocation = "{ICO},0"
$s.Description = "JUICE - chamber experiment to model run"
$s.Save()
'''
    r = subprocess.run(["powershell", "-NoProfile", "-Command", ps], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"could not create {link_path}: {r.stderr.strip()}")
    print("created", link_path)
    set_app_id(link_path)


APP_ID = "FZJ.JUICE"      # the AppUserModelID every JUICE window sets (juice_theme.taskbar_identity)


def set_app_id(link_path: str, app_id: str = APP_ID) -> bool:
    """Stamp the AppUserModelID into the shortcut: Windows then matches the running JUICE windows
    to this shortcut and shows its JUICE icon on the taskbar (pinned or not) instead of python's.
    Needs pywin32; silently skipped without it."""
    try:
        import pythoncom
        from win32com.propsys import propsys, pscon
        GPS_READWRITE = 2
        ps = propsys.SHGetPropertyStoreFromParsingName(link_path, None, GPS_READWRITE, propsys.IID_IPropertyStore)
        ps.SetValue(pscon.PKEY_AppUserModel_ID, propsys.PROPVARIANTType(app_id, pythoncom.VT_LPWSTR))
        ps.Commit()
        print("  AppUserModelID", app_id, "written")
        return True
    except Exception as e:
        print("  (AppUserModelID not written:", e, ")")
        return False


if __name__ == "__main__":
    if not sys.platform.startswith("win"):
        raise SystemExit("Windows only")
    if "--here" in sys.argv:
        make(os.path.join(HERE, "JUICE.lnk"))
    else:
        desktop = os.path.join(os.path.expanduser("~"), "Desktop")
        make(os.path.join(desktop, "JUICE.lnk"))
    if "--startmenu" in sys.argv:
        start = os.path.join(os.environ.get("APPDATA", ""), "Microsoft", "Windows", "Start Menu", "Programs")
        make(os.path.join(start, "JUICE.lnk"))
