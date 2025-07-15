import streamlit as st
from semantic_search import search_documents, generate_summary

st.set_page_config(page_title="🧠 PubMed Semantic Search", layout="wide")
st.title("🔍 Medical Journals Semantic Search")

# Input fields
query = st.text_input("📝 Enter your medical query:")
top_k = st.number_input("📄 Number of relevant documents (top_k):", min_value=1, max_value=20, value=5)

# Submit
if st.button("🔍 Search PubMed Abstracts"):
    if not query.strip():
        st.warning("Please enter a valid query.")
    else:
        with st.spinner("Processing..."):
            results = search_documents(query, top_k=top_k)

            st.subheader("📚 Top Retrieved Documents")
            for i, doc in enumerate(results, 1):
                meta = doc.metadata or {}
                title = meta.get("title", "N/A")
                pmid = meta.get("pmid", "N/A")
                authors = ", ".join(meta.get("authors", [])) if meta.get("authors") else "N/A"
                journal = meta.get("journal", "N/A")
                url = meta.get("url", "")
                score = doc.metadata.get("score", None)

                st.markdown(f"### 📄 Document {i}")
                st.markdown(f"**Score:** {score:.4f}" if score else "Score: N/A")
                st.markdown(f"**Title:** {title}")
                st.markdown(f"**PMID:** {pmid}")
                st.markdown(f"**Authors:** {authors}")
                st.markdown(f"**Journal:** {journal}")
                if url:
                    st.markdown(f"🔗 [View on PubMed]({url})")
                st.markdown(f"> {doc.page_content or 'No content available.'}")
                st.markdown("---")

            st.subheader("🧠 Gemini Summary Answer")
            answer = generate_summary(query, results)
            st.success(answer)
