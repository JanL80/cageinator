
import sys, json



def dumps_inline(obj, indent=2, inline_keys={"xyz", "vector"}):
    def fmt(v, level=0, parent_key=None):
        sp = " " * (indent * level)
        if isinstance(v, dict):
            if not v: return "{}"
            items = list(v.items())
            lines = ["{"]
            for i, (k, vv) in enumerate(items):
                comma = "," if i < len(items) - 1 else ""
                lines.append(
                    " " * (indent * (level + 1))
                    + json.dumps(k) + ": "
                    + fmt(vv, level + 1, parent_key=k)
                    + comma
                )
            lines.append(sp + "}")
            return "\n".join(lines)
        if isinstance(v, list):
            if parent_key in inline_keys and all(not isinstance(x, (list, dict)) for x in v):
                return "[" + ", ".join(json.dumps(x) for x in v) + "]"
            if not v: return "[]"
            lines = ["["]
            for i, x in enumerate(v):
                comma = "," if i < len(v) - 1 else ""
                lines.append(" " * (indent * (level + 1)) + fmt(x, level + 1) + comma)
            lines.append(sp + "]")
            return "\n".join(lines)
        return json.dumps(v)
    return fmt(obj, 0, None) + "\n"



def pick_name():

    if len(sys.argv) > 1 and sys.argv[1] != "--":           ### CLI arg
        return sys.argv[1]

    try:                                                    ### piped stdin (non-blocking when TTY)
        if not sys.stdin.isatty():
            data = sys.stdin.read().strip()
            if data:
                return data
    except Exception:
        pass

    return "temp"                                           ### default



def main():
    name = pick_name()
    if not name.endswith(".json"):
        name += ".json"

    skeleton = {
      "unit_id": "",
      "class": "",
      "gbu_type": "",
      "gbu_subtype": "",
      "composition": {"formula": "", "charge": None, "spin_mult": None},
      "atoms": [{"idx": None, "el": "", "xyz": [None, None, None]}],
      "bonds": [{"a": None, "b": None, "order": None}],
      "connectors": [{
        "id": "",
        "atom_index": None,
        "element": "",
        "vector": [None, None, None],
        "role": "",
        "site_label": "",
        "constraints": {"preferred_partner": ""}
      }],
      "coordination_count": None,
      "coordination_atoms": [None]
    }

    with open(name, "w") as f:
        f.write(dumps_inline(skeleton, inline_keys={"xyz", "vector"}))
    print(name)



if __name__ == "__main__":
    main()
