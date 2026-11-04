r"""
make_shortcut.py — put a JUICE icon on the Desktop (and optionally the Start menu).

    python make_shortcut.py              # Desktop\JUICE.lnk → JUICE.bat with the JUICE icon
    python make_shortcut.py --here       # JUICE.lnk inside the project folder instead
    python make_shortcut.py --startmenu  # also Start Menu\Programs\JUICE.lnk

One icon for everything: it opens the JUICE start menu (juice_start.py): new experiment, saved
experiment, new study, results of one experiment, results of a study.  Each click opens a new,
independent menu (JUICE.bat uses `start` + pythonw), so another window can be opened while others
are running.  A file dropped on the icon opens straight in the right window.  The separate
"JUICE Results" icon of earlier versions is removed when this script finds it.  Windows only.

A .lnk always stores absolute paths of the machine it was made on, so it is NOT synced
through git (JUICE.lnk is ignored); JUICE.bat is, and works anywhere.  After a clone or on
another PC run this script once to get the icon.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BAT = os.path.join(HERE, "JUICE.bat")
BAT_RESULTS = os.path.join(HERE, "JUICE Results.bat")
ICO = os.path.join(HERE, "JUICE.ico")


DESCRIPTION = "JUICE - set up and run chamber experiments and studies, or look at their results"


def remove_old_results_icon(folder: str) -> None:
    """Delete the 'JUICE Results' shortcut earlier versions of this script made (only if it points to our .bat)."""
    old = os.path.join(folder, "JUICE Results.lnk")
    if not os.path.isfile(old):
        return
    ps = f'(New-Object -ComObject WScript.Shell).CreateShortcut("{old}").TargetPath'
    r = subprocess.run(["powershell", "-NoProfile", "-Command", ps], capture_output=True, text=True)
    if os.path.normcase(r.stdout.strip()) == os.path.normcase(BAT_RESULTS):
        os.remove(old)
        print("removed", old, "(the start menu of the JUICE icon covers it)")


def make(link_path: str, target: str = BAT, description: str = DESCRIPTION) -> None:
    ps = f'''
$s = (New-Object -ComObject WScript.Shell).CreateShortcut("{link_path}")
$s.TargetPath = "{target}"
$s.WorkingDirectory = "{HERE}"
$s.IconLocation = "{ICO},0"
$s.Description = "{description}"
$s.Save()
'''
    r = subprocess.run(["powershell", "-NoProfile", "-Command", ps], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"could not create {link_path}: {r.stderr.strip()}")
    print("created", link_path)
    set_app_id(link_path)


def desktop_folder() -> str:
    """The real Desktop (OneDrive can move it away from ~/Desktop), asked from Windows directly
    (a PowerShell round trip garbles non-ASCII paths such as 'Jülich')."""
    try:
        import ctypes
        buf = ctypes.create_unicode_buffer(32768)
        if ctypes.windll.shell32.SHGetFolderPathW(None, 0x10, None, 0, buf) == 0 and os.path.isdir(buf.value):  # CSIDL_DESKTOPDIRECTORY
            return buf.value
    except Exception:
        pass
    return os.path.join(os.path.expanduser("~"), "Desktop")


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
    folders = [HERE if "--here" in sys.argv else desktop_folder()]
    if "--startmenu" in sys.argv:
        folders.append(os.path.join(os.environ.get("APPDATA", ""), "Microsoft", "Windows", "Start Menu", "Programs"))
    for f in folders:
        make(os.path.join(f, "JUICE.lnk"))
        remove_old_results_icon(f)
