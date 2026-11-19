#!/usr/bin/env python3
"""Download CTRPv2 etc. drug-response profiles from PharmacoDB."""
import json, os, time, urllib.request
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cache", "pharmacodb")
URL = "https://pharmacodb.ca/graphql"
def gql(query):
    req = urllib.request.Request(URL, data=json.dumps({"query": query}).encode(), headers={"Content-Type": "application/json", "User-Agent": "python-urllib"})
    err = None
    for attempt in range(3):
        try:
            return json.load(urllib.request.urlopen(req, timeout=240))
        except Exception as ex:
            err = str(ex); time.sleep(3)
    return {"errors": err}
UIDS = {"ML210": "PDBC00340", "ML162": "PDBC00312", "RSL3": "PDBC00401", "Erastin": "PDBC00307", "Auranofin": "PDBC31429",
        "BSO": "PDBC02790", "Sulfasalazine": "PDBC00601", "Piperlongumine": "PDBC00402"}
for name, uid in UIDS.items():
    r = gql('{ compound(compoundUID: "%s") { compound { id name } } }' % uid)
    try: cid = r["data"]["compound"]["compound"]["id"]
    except Exception: print(name, "no id", str(r)[:200]); continue
    r = gql('{ experiments(compoundId: %d, all: true) { cell_line { name } tissue { name } dataset { name } profile { AAC IC50 } } }' % cid)
    ex = (r.get("data") or {}).get("experiments")
    if ex is None:
        print(name, cid, "ERROR", str(r)[:300], flush=True); continue
    json.dump(ex, open(f"{OUT}/experiments_{name}.json", "w"))
    ds = {}
    for e in ex: ds[e["dataset"]["name"]] = ds.get(e["dataset"]["name"], 0) + 1
    print(name, cid, len(ex), "experiments;", ds, flush=True)
