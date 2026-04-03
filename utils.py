import itertools as it

def _coerce(value):
    """Cast a string override value to the appropriate Python type."""
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    if value.lower() in ("null", "none"):
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


def apply_overrides(args, overrides):
    """
    Apply a list of 'key=value' or 'nested.key=value' override strings to an
    args dict, coercing types automatically.

    Examples
    --------
    --set iter=5
    --set pedigree.mating=mono
    --set run_ibdne=false
    --set custom_demo.object=ooa2
    """
    for item in it.chain(*overrides):
        if "=" not in item:
            raise ValueError(
                f"Invalid override '{item}' — expected format KEY=VALUE or NESTED.KEY=VALUE"
            )
        raw_key, raw_val = item.split("=", 1)
        keys = raw_key.strip().split(".")
        value = _coerce(raw_val.strip())

        # Walk / create nested dicts
        d = args
        for k in keys[:-1]:
            if k not in d or not isinstance(d[k], dict):
                d[k] = {}
            d = d[k]
        d[keys[-1]] = value
        print(f"Override: {raw_key} = {value!r}")

    return args
