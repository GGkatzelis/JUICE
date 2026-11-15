#!/usr/bin/env python
"""parse_cache.py -- the parse cache of scenario_runner (ScenarioRunner._parse_elog_cached), moved out
unchanged 2026-11-15.

    m_part = parse_cached(parser, elog_path, enabled=True, log=print, source_file=None)

parser.parse(elog_path), with the parser state of large elog files (>= 5 MB) kept in ~/.juice/parse_cache.
The key is a hash of the elog text (includes resolved), the parser sources (easy_parser.py, mechanism.py,
evaluator.py in the directory of `source_file`) and the Python version, so any change of the inputs or of the
parser gives a new entry.  The last six entries are kept.

`source_file`: a file in the directory that holds the three parser sources.  scenario_runner passes its own
__file__, so the key is the one the method computed before the move; None = the directory of this file.
"""
import os
import sys
import time

from juice.core.mechanism import Mechanism


def parse_cached(parser, elog_path, enabled=True, log=print, source_file=None):
    """parser.parse(elog_path), with the parser state of large elog files kept in ~/.juice/parse_cache.
    The key is a hash of the elog text (includes resolved), the parser sources and the Python version, so
    any change of the inputs or of the parser gives a new entry.  `enabled = False` switches it off."""
    _on = enabled
    _f = None
    try:
        if _on and os.path.getsize(elog_path) >= 5_000_000:
            import hashlib, pickle
            _t0 = time.time()
            _content = parser.preprocess(elog_path)
            _h = hashlib.sha256(_content.encode('utf-8', 'surrogatepass'))
            _here = os.path.dirname(os.path.abspath(source_file or __file__))
            for _m in ('easy_parser.py', 'mechanism.py', 'evaluator.py'):
                with open(os.path.join(_here, _m), 'rb') as _fh:
                    _h.update(_fh.read())
            _h.update(sys.version.encode()); _h.update(b'elog-cache-v2')
            _dir = os.path.join(os.path.expanduser('~'), '.juice', 'parse_cache')
            os.makedirs(_dir, exist_ok=True)
            _f = os.path.join(_dir, _h.hexdigest()[:32] + '.pkl')
            if os.path.isfile(_f):
                with open(_f, 'rb') as _fh:
                    _state = pickle.load(_fh)
                _plog = _state.pop('__parse_log__', '')
                if not isinstance(_state.get('mech'), Mechanism):
                    raise ValueError('cache entry without mechanism')
                try:
                    os.utime(_f, None)
                except Exception:
                    pass
                try:
                    log(f"  [MECHANISM] parsed mechanism taken from the parse cache ({time.time() - _t0:.1f} s; "
                        f"{os.path.basename(_f)})")
                    if _plog.strip():
                        log("  [MECHANISM] messages of the parser when this entry was written:")
                        log(_plog.rstrip())
                except Exception:
                    pass
                # last step, nothing after it can fail: a failure above falls back to a fresh parse of the
                # still empty parser
                parser.__dict__.update(_state)
                return parser.mech
            parser.preprocess = lambda _fn, _c=_content: _c        # parse() reads the text once more: reuse it
    except Exception as _e:
        log(f"  [MECHANISM] parse cache not used ({_e})")
        _f = None
        parser.__dict__.pop('preprocess', None)
    import io as _io, contextlib as _ctxlib
    _buf = _io.StringIO()
    try:
        if _f:      # keep the parser's messages with the cache entry (they are printed again on a hit)
            with _ctxlib.redirect_stdout(_buf):
                m_part = parser.parse(elog_path)
        else:
            m_part = parser.parse(elog_path)
    finally:
        parser.__dict__.pop('preprocess', None)
        if _buf.getvalue():
            log(_buf.getvalue(), end='')
    if _f:
        try:
            import pickle
            _state = {k: v for k, v in parser.__dict__.items() if k not in ('base_path', 'include_dir', 'preprocess')}
            _state['__parse_log__'] = _buf.getvalue()
            _tmp = f"{_f}.{os.getpid()}.tmp"
            with open(_tmp, 'wb') as _fh:
                pickle.dump(_state, _fh, protocol=pickle.HIGHEST_PROTOCOL)
            os.replace(_tmp, _f)
            _dir = os.path.dirname(_f)
            _old = sorted((os.path.join(_dir, x) for x in os.listdir(_dir) if x.endswith('.pkl')), key=os.path.getmtime)
            for _x in _old[:-6]:
                os.remove(_x)
            log(f"  [MECHANISM] parse cache written ({os.path.getsize(_f) / 1e6:.0f} MB, {_dir})")
        except Exception as _e:
            log(f"  [MECHANISM] parse cache not written ({_e})")
    return m_part
