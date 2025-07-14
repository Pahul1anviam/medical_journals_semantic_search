# import os
# import numpy as np
# from dotenv import load_dotenv
# import google.generativeai as genai
# from sentence_transformers import SentenceTransformer
# from qdrant_client import QdrantClient

# # Load Gemini key
# load_dotenv()
# genai.configure(api_key=os.getenv("GOOGLE_API_KEY"))

# # Init models
# embedder = SentenceTransformer("BAAI/bge-small-en-v1.5")
# client = QdrantClient(path="./qdrant_db")
# collection_name = "medical_chunks"

# # Embed query
# def embed(text):
#     return embedder.encode([text], convert_to_numpy=True)[0]

# # Semantic search
# def semantic_search(query, top_k=5):
#     query_vec = embed(query)
#     results = client.search(
#         collection_name=collection_name,
#         query_vector=query_vec,
#         limit=top_k
#     )
#     return results

# # Ask Gemini to summarize
# def generate_answer(query, results):
#     context = "\n\n".join([r.payload["text"] for r in results])
#     prompt = f"""
# You are a clinical research assistant.

# Answer this question: **{query}**
# Use only the context below, no external information.

# Context:
# {context}
# """
#     model = genai.GenerativeModel("gemini-2.5-flash")
#     response = model.generate_content(prompt)
#     return response.text

# # --- Run ---
# if __name__ == "__main__":
#     query = "What are the treatment options for pancreatic cancer in elderly patients?"
#     results = semantic_search(query)

#     print("🔍 Top Results:")
#     for r in results:
#         print(f"\nScore: {r.score:.4f}")
#         print(r.payload["text"][:300])
#         print("-" * 50)

#     print("\n🧠 Gemini Response:\n")
#     print(generate_answer(query, results))









# semantic_search.py
import os
import numpy as np
from dotenv import load_dotenv
import google.generativeai as genai
from sentence_transformers import SentenceTransformer
from qdrant_client import QdrantClient

# Load environment and Gemini API key
load_dotenv()
genai.configure(api_key=os.getenv("GOOGLE_API_KEY"))

# Initialize models and client
embedder = SentenceTransformer("BAAI/bge-small-en-v1.5")
client = QdrantClient(path="./qdrant_db")
collection_name = "medical_chunks"

# Embed text
def embed(text):
    return embedder.encode([text], convert_to_numpy=True)[0]

# Semantic search function
def semantic_search(query, top_k=5):
    query_vec = embed(query)
    results = client.search(
        collection_name=collection_name,
        query_vector=query_vec,
        limit=top_k
    )
    return results

# Gemini: Format each individual chunk
def summarize_chunk(text):
    prompt = f"""
Extract the key findings from the following medical research snippet and format them as clear bullet points:

{text}
"""
    model = genai.GenerativeModel("gemini-2.5-flash")
    response = model.generate_content(prompt)
    return response.text

# Gemini: Summarize all results together in structured format
def generate_summary(query, results):
    context = "\n\n".join([r.payload["text"] for r in results])
    prompt = f"""
You are a clinical research assistant.

Using only the information provided below, summarize the answer to the question in the following structure:

- ✅ Key Findings:
- 📊 Notable Statistics:
- 💡 Clinical Implications:
- 🧪 Limitations or Future Research Suggestions:

Question:
{query}

Context:
{context}
"""
    model = genai.GenerativeModel("gemini-2.5-flash")
    response = model.generate_content(prompt)
    return response.text

# --- Run as script ---
if __name__ == "__main__":
    query = "What are the common complications of diabetes in elderly patients?"
    top_k = 5

    print(f"🔍 Searching for: {query}")
    results = semantic_search(query, top_k)

    print("\n📚 Top Retrieved Chunks:\n")
    for i, r in enumerate(results, 1):
        print(f"\n{i}. Score: {r.score:.4f}")
        print(summarize_chunk(r.payload["text"]))
        print("-" * 80)

    print("\n🧠 Gemini Final Structured Summary:\n")
    print(generate_summary(query, results))
