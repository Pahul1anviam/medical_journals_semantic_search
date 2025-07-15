from Bio import Entrez, Medline
from tqdm import tqdm
import time
import json
import ssl

ssl._create_default_https_context = ssl._create_unverified_context
Entrez.email = "pahuljot.singh@anviam.com"

def fetch_pubmed_ids(query, max_results=10000):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    return record["IdList"]

def fetch_abstracts_by_ids(id_list):
    abstracts = []
    for i in tqdm(range(0, len(id_list), 100), desc="Fetching abstracts"):
        chunk = id_list[i:i+100]
        ids = ",".join(chunk)
        try:
            handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
            records = Medline.parse(handle)
            for record in records:
                abstracts.append({
                    "pmid": record.get("PMID", ""),
                    "title": record.get("TI", ""),
                    "authors": record.get("AU", []),
                    "abstract": record.get("AB", ""),
                    "journal": record.get("JT", ""),
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{record.get('PMID', '')}/"
                })
            time.sleep(0.3)
        except Exception as e:
            print(f"Error: {e}")
            continue
    return abstracts

def save_to_jsonl(data, output_file="pubmed_abstracts_new.jsonl"):
    with open(output_file, "w", encoding="utf-8") as f:
        for item in data:
            if len(item.get("abstract", "").strip()) > 100:
                f.write(json.dumps(item) + "\n")

if __name__ == "__main__":
    query = "diabetes OR cancer OR covid OR heart disease OR lung disease"
    ids = fetch_pubmed_ids(query)
    print(f"Found {len(ids)} IDs")
    data = fetch_abstracts_by_ids(ids)
    save_to_jsonl(data)
    print(f"Saved to pubmed_abstracts.jsonl")
