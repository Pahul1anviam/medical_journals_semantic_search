import os
import json
import numpy as np
from tqdm import tqdm
from uuid import uuid4
from dotenv import load_dotenv
from sentence_transformers import SentenceTransformer
from qdrant_client import QdrantClient
from qdrant_client.models import PointStruct


load_dotenv()

# Initialize
client = QdrantClient(path="./qdrant_db")
collection_name = "medical_chunks"
model = SentenceTransformer("BAAI/bge-small-en-v1.5")

# Load chunks
chunks = []
with open("pubmed_chunks.jsonl", "r", encoding="utf-8") as f:
    for line in f:
        obj = json.loads(line)
        chunks.append(obj["chunk"])

print(f"✅ Loaded {len(chunks)} chunks")

# Embed and insert
def embed(texts):
    if isinstance(texts, str):
        texts = [texts]
    return model.encode(texts, convert_to_numpy=True)

batch_size = 50
for i in tqdm(range(0, len(chunks), batch_size), desc="Embedding + Uploading"):
    batch = chunks[i:i+batch_size]
    embeddings = embed(batch)

    client.upsert(
    collection_name=collection_name,
    points=[
        PointStruct(
            id=str(uuid4()),
            vector=emb.tolist(),
            payload={
                "text": batch[j],
                "source": f"paper-{i+j}"
            }
        )
        for j, emb in enumerate(embeddings)
    ]
)

print("✅ All chunks embedded and inserted")
