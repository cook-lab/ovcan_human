#!/usr/bin/env python3
"""
build_single.py — produce a SINGLE self-contained HTML for deployment.

Inlines app/data.js into app/index.html (replacing the external <script src>),
so the result has zero external references: it renders identically at
cooklab.ca/ovcan_viewer, /ovcan_viewer/, or from file:// — no trailing-slash /
relative-path fragility, no server config.

Usage:  python3 app/build_single.py
Output: app/ovcan_viewer_standalone.html  (copy to the site as .../index.html)

The two-file version (index.html + data.js) stays the source of truth for local
dev; regenerate data.js with build_payload.py, then re-run this to refresh the
standalone. base64 payload contains no '<', so inlining cannot break the tag.
"""
import os

APP = os.path.dirname(os.path.abspath(__file__))
html = open(os.path.join(APP, "index.html"), encoding="utf-8").read()
data = open(os.path.join(APP, "data.js"), encoding="utf-8").read()

marker = '<script src="data.js"></script>'
if marker not in html:
    raise SystemExit("ERROR: expected external include %r not found in index.html" % marker)
if "</script" in data.lower():
    raise SystemExit("ERROR: data.js contains a </script> sequence; inlining would break the tag")

html = html.replace(marker, "<script>\n" + data + "\n</script>")

out = os.path.join(APP, "ovcan_viewer_standalone.html")
open(out, "w", encoding="utf-8").write(html)
print("wrote %s  (%.2f MB, self-contained)" % (out, len(html.encode("utf-8")) / 1e6))
