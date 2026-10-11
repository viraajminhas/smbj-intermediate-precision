"""Download the localisation tables used by analyse.py from BioImage Archive S-BIAD8: the U2OS
Nup96-SNAP-AF647 cells imaged in the Gloxy + 35 mM MEA buffer. Reads the file listing smlm_data.tsv;
saves to data/ and skips files already present.
"""
import urllib.parse, urllib.request, os, csv, re, collections
base="https://www.ebi.ac.uk/biostudies/files/S-BIAD8/"
rows=list(csv.DictReader(open('smlm_data.tsv',encoding='utf-8'),delimiter='\t'))
pat=re.compile(r'_(\d{6})_')
# matched group: U2OS Nup96 SNAP-AF647, identical buffer folder, multiple dates
want=[r['Files'] for r in rows if r['Files'].endswith('_sml.csv')
      and 'SNAP/BufferTesting/Gloxy+35 mM MEA' in r['Files']]
print("files in matched buffer group:",len(want))
os.makedirs("data",exist_ok=True)
for f in want:
    out="data/"+os.path.basename(f)
    if os.path.exists(out) and os.path.getsize(out)>1000:
        print("  have",os.path.basename(f)); continue
    url=base+urllib.parse.quote(f)
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'})
        d=urllib.request.urlopen(req,timeout=300).read()
        open(out,"wb").write(d)
        print(f"  got {os.path.basename(f)} {len(d)/1e6:.1f} MB")
    except Exception as e:
        print("  ERR",os.path.basename(f),e)
