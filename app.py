import streamlit as st
from qdrant_client import QdrantClient
from sentence_transformers import SentenceTransformer
from langchain.vectorstores import Qdrant
from langchain.embeddings import HuggingFaceEmbeddings
from langchain.llms import GoogleGenerativeAI
from langchain.chains import RetrievalQA
from langchain.prompts import PromptTemplate
import os

# ----------------------------- Config ----------------------------- #
QDRANT_URL = "http://localhost:6333"
COLLECTION_NAME = "pubmed_chunks"

embedding_model = HuggingFaceEmbeddings(model_name="sentence-transformers/all-MiniLM-L6-v2")
llm = GoogleGenerativeAI(model="gemini-1.5-flash", google_api_key=os.getenv("GOOGLE_API_KEY"))

# -------------------------- Streamlit UI -------------------------- #
st.set_page_config(page_title="PubMed Semantic Search", layout="wide")
st.title("🔍 PubMed Semantic Search Assistant")

st.markdown("""
Enter a medical research query to retrieve the most relevant PubMed abstracts and get an AI-generated summary.
""")

# Sample Query
st.caption("💡 *Example: What are recent findings on mRNA vaccine effectiveness in elderly populations?*")

query = st.text_input("Enter your research query:", value="")

top_k = st.slider("Number of relevant documents to retrieve", min_value=1, max_value=10, value=3)

if st.button("Search PubMed"):
    if not query:
        st.warning("Please enter a query.")
    else:
        # ---------------------- Load Vector Store ---------------------- #
        qdrant = Qdrant(
            client=QdrantClient(url=QDRANT_URL),
            collection_name=COLLECTION_NAME,
            embeddings=embedding_model,
        )

        retriever = qdrant.as_retriever(search_kwargs={"k": top_k})

        # --------------------- Custom Prompt --------------------- #
        prompt_template = PromptTemplate(
            input_variables=["context", "question"],
            template="""
You are a biomedical research assistant. Based on the provided context from medical research articles, generate a concise summary.

Context:
{context}

Question:
{question}

Answer in clinical terms:
""")

        chain = RetrievalQA.from_chain_type(
            llm=llm,
            retriever=retriever,
            return_source_documents=True,
            chain_type_kwargs={"prompt": prompt_template},
        )

        # ---------------------- Run Chain ---------------------- #
        response = chain(query)

        st.subheader("🧠 Gemini Summary")
        st.markdown(response["result"])

        # ---------------------- Show Documents ---------------------- #
        with st.expander("📄 Show Retrieved Chunks"):
            for i, doc in enumerate(response["source_documents"]):
                st.markdown(f"**Chunk {i+1}**")
                st.code(doc.page_content, language="markdown")
