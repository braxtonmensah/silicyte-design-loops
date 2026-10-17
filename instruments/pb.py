"""Porkbun API helper. Reads credentials from a file and never prints them."""
import json, sys, urllib.request

CRED = json.load(open(r"C:/Users/bmens/.porkbun_api.json"))
BASE = "https://api.porkbun.com/api/json/v3/"

def call(path, **kw):
    body = dict(CRED); body.update(kw)
    req = urllib.request.Request(BASE + path,
        data=json.dumps(body).encode(), headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=45) as r:
            return json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        body = e.read().decode(errors="replace")
        try:
            return json.loads(body)
        except Exception:
            return {"status": "HTTPERROR", "code": e.code, "body": body[:500]}

if __name__ == "__main__":
    print(json.dumps(call(sys.argv[1], **dict(a.split("=",1) for a in sys.argv[2:])), indent=1)[:4000])
