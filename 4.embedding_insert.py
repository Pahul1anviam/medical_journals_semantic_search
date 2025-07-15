import json
from uuid import uuid4
from tqdm import tqdm
import numpy as np
from qdrant_client import QdrantClient
from qdrant_client.models import PointStruct
from sentence_transformers import SentenceTransformer

client = QdrantClient(path="./qdrant_db")
model = SentenceTransformer("BAAI/bge-small-en-v1.5")

with open("pubmed_chunks_new.jsonl", "r", encoding="utf-8") as f:
    records = [json.loads(line) for line in f]

print(f"✅ Loaded {len(records)} records")

def embed(texts):
    return model.encode(texts, convert_to_numpy=True)

for i in tqdm(range(0, len(records), 50), desc="Embedding and Uploading"):
    batch = records[i:i+50]
    texts = [r["chunk"] for r in batch]
    embeddings = embed(texts)

    points = []
    for j, emb in enumerate(embeddings):
        record = batch[j]
        points.append(PointStruct(
            id=str(uuid4()),
            vector=emb.tolist(),
            payload={
                "text": record["chunk"],
                "pmid": record["pmid"],
                "title": record["title"],
                "authors": record["authors"],
                "url": record["url"],
                "journal": record["journal"]
            }
        ))
    client.upsert(collection_name="medical_chunks_new", points=points)

print("✅ Done")
