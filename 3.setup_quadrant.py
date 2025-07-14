from qdrant_client import QdrantClient
from qdrant_client.models import VectorParams, Distance

collection_name = "medical_chunks"

client = QdrantClient(path="./qdrant_db")
client.recreate_collection(
    collection_name=collection_name,
    vectors_config=VectorParams(size=384, distance=Distance.COSINE)
)

print("✅ Qdrant collection created")
