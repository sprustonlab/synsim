"""Emit the simulation parameter spec as JSON for the front-end to build controls.

    python scripts/dump_params.py            # print to stdout
    python scripts/dump_params.py web/params.json   # write to a file
"""
import json
import sys

from synsim import ui_spec, DEFAULTS


def main():
    doc = {"defaults": DEFAULTS.to_dict(), "controls": ui_spec()}
    text = json.dumps(doc, indent=2)
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w", encoding="utf-8") as f:
            f.write(text)
        print(f"-> {sys.argv[1]} ({len(doc['controls'])} controls)")
    else:
        print(text)


if __name__ == "__main__":
    main()
