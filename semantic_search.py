import os
from dotenv import load_dotenv
from qdrant_client import QdrantClient
from qdrant_client.models import Filter, FieldCondition, MatchValue
from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain.prompts import PromptTemplate
from langchain.chains import LLMChain
from langchain.schema import Document

# Load environment variables and API key
load_dotenv()

# Initialize embedding model
embedder = HuggingFaceEmbeddings(model_name="BAAI/bge-small-en-v1.5")

# Initialize Qdrant client
client = QdrantClient(path="./qdrant_db")

# Initialize Gemini model
llm = ChatGoogleGenerativeAI(model="gemini-2.5-flash")

# Prompt template for summarization
SUMMARY_PROMPT = PromptTemplate(
    input_variables=["query", "context"],
    template="""
You are a clinical research assistant.

Using only the provided content, answer the question below using this structure:

- ✅ Key Findings
- 📊 Notable Statistics
- 💡 Clinical Implications
- 🔗 Metadata (title, authors, journal, pmid, and url)

Question:
{query}

Context:
{context}
"""
)

# Gemini Chain
chain = LLMChain(llm=llm, prompt=SUMMARY_PROMPT)

# ----------- 🔍 Semantic Search Using Raw Qdrant --------------
def search_documents(query, top_k=5):
    query_vector = embedder.embed_query(query)
    results = client.search(
        collection_name="medical_chunks_new",
        query_vector=query_vector,
        limit=top_k,
        with_payload=True,
        with_vectors=False
    )

    docs = []
    for r in results:
        doc = Document(
            page_content=r.payload.get("text", ""),
            metadata={**r.payload, "score": r.score}  # ⬅ Store score inside metadata
        )
        docs.append(doc)
    return docs

# ----------- 🧠 Gemini Summary Generation ---------------------
def generate_summary(query, docs):
    context = "\n\n".join(doc.page_content for doc in docs)
    return chain.run({"query": query, "context": context})
