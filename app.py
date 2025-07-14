# app.py

import streamlit as st
from semantic_search import semantic_search, generate_summary


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
            results = semantic_search(query, top_k=top_k)
            st.subheader("📚 Top Retrieved Chunks")
            for i, r in enumerate(results, 1):
                st.markdown(f"**{i}. Score:** {r.score:.4f}")
                st.write(r.payload["text"])
                st.markdown("---")

            st.subheader("🧠 Gemini Summary Answer")
            answer = generate_summary(query, results)
            st.success(answer)
