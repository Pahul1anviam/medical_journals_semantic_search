import streamlit as st
from semantic_search import search_documents, generate_summary

st.set_page_config(page_title="Medical Research Navigator", layout="wide")
st.title("🩺 Medical Research Navigator")
st.markdown("Semantic search + Gemini-based summarization for PubMed abstracts.")

# 🔢 Select how many documents (chunks) to retrieve
top_k = st.slider("📄 Number of relevant documents to fetch", min_value=1, max_value=10, value=5)

# 🔍 User query
query = st.text_input("Enter your clinical query", placeholder="e.g., effectiveness of metformin in Indian patients")

if st.button("Search"):
    if query.strip() == "":
        st.warning("Please enter a valid query.")
    else:
        with st.spinner("🔍 Searching and summarizing..."):
            docs = search_documents(query, top_k=top_k)  
            answer = generate_summary(query, docs)

        st.markdown("### Gemini Summary")
        st.markdown(answer)

        st.markdown("---")
        st.markdown(f"### 📄 Top {len(docs)} Source Chunks Used")

        for i, doc in enumerate(docs):
            meta = doc.metadata
            st.markdown(f"#### 📄 Document {i+1}")
            st.markdown(f"**Score:** {meta.get('score', 0):.4f}")
            st.markdown(f"**Title:** {meta.get('title', 'N/A')}")
            st.markdown(f"**Journal:** {meta.get('journal', 'N/A')}")
            st.markdown(f"**Authors:** {', '.join(meta.get('authors', []))}")
            st.markdown(f"**PMID:** {meta.get('pmid', 'N/A')}")
            st.markdown(f"🔗 [View on PubMed]({meta.get('url', '#')})")

            with st.expander("Show Source Chunk"):
                st.markdown(doc.page_content)

            st.markdown("---")
