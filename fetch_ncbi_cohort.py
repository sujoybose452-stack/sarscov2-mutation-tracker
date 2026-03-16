import argparse
import time
import urllib.parse
import urllib.request
import ssl
import xml.etree.ElementTree as ET


def _open_url(url, insecure=False):
    if insecure:
        ctx = ssl._create_unverified_context()
        return urllib.request.urlopen(url, context=ctx)
    return urllib.request.urlopen(url)


def esearch(term, retmax=100, db="nuccore", insecure=False):
    params = {
        "db": db,
        "term": term,
        "retmax": str(retmax),
        "retmode": "xml",
        "tool": "codex",
    }
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
    with _open_url(url, insecure=insecure) as resp:
        xml_text = resp.read()
    root = ET.fromstring(xml_text)
    ids = [elem.text for elem in root.findall(".//IdList/Id")]
    return ids


def efetch(ids, db="nuccore", rettype="fasta", insecure=False):
    params = {
        "db": db,
        "id": ",".join(ids),
        "rettype": rettype,
        "retmode": "text",
        "tool": "codex",
    }
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
    with _open_url(url, insecure=insecure) as resp:
        return resp.read().decode("utf-8")


def main():
    parser = argparse.ArgumentParser(description="Fetch SARS-CoV-2 cohort FASTA from NCBI")
    parser.add_argument("--query", required=True, help="NCBI search query")
    parser.add_argument("--retmax", type=int, default=50, help="Number of sequences to fetch")
    parser.add_argument("--out", required=True, help="Output FASTA path")
    parser.add_argument("--chunk", type=int, default=200, help="IDs per efetch request")
    parser.add_argument("--sleep", type=float, default=0.34, help="Seconds between efetch calls")
    parser.add_argument(
        "--insecure",
        action="store_true",
        help="Disable SSL verification if your environment blocks HTTPS",
    )
    args = parser.parse_args()

    ids = esearch(args.query, retmax=args.retmax, insecure=args.insecure)
    if not ids:
        raise SystemExit("No IDs returned from NCBI search. Try adjusting the query.")

    chunks = [ids[i : i + args.chunk] for i in range(0, len(ids), args.chunk)]
    with open(args.out, "w") as handle:
        for i, chunk in enumerate(chunks, start=1):
            fasta = efetch(chunk, insecure=args.insecure)
            handle.write(fasta)
            if i < len(chunks):
                time.sleep(args.sleep)


if __name__ == "__main__":
    main()
