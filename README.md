Medical Journals Semantic Search
================================

This project is a semantic search engine for medical research papers, built using local embeddings and Google Gemini LLM.

It allows users to ask medical questions in plain English and retrieves the most relevant research papers based on semantic similarity. A final summarized answer is also generated using Gemini.

-------------------------------------
Features
--------

- Search medical papers using natural language
- Uses Hugging Face embeddings (bge-small-en-v1.5)
- Stores vectors locally using Qdrant (no Docker)
- Gemini generates answer summaries from results
- Streamlit dashboard to run the app

-------------------------------------
Project Structure
-----------------

1.fetch.py              - Downloads PubMed-style abstracts  
2.tokenization.py       - Splits abstracts into chunks  
3.setup_qdrant.py       - Sets up the local Qdrant DB  
4.embedding_insert.py   - Embeds chunks and stores in DB  
5.semantic_search.py    - Runs semantic search + LLM summary (CLI)  
app.py                  - Streamlit UI  
pubmed_abstracts.jsonl  - Raw abstract data  
pubmed_chunks.jsonl     - Preprocessed chunks  
qdrant_db/              - Local Qdrant vector storage  
.env                    - Store your Google Gemini API key  
requirements.txt        - Python dependencies

-------------------------------------
How to Run
----------

1. Create a virtual environment:
   > python -m venv venv
   > venv\Scripts\activate   (for Windows)

2. Install dependencies:
   > pip install -r requirements.txt

3. Add your Gemini API key in a `.env` file:
   GOOGLE_API_KEY=your_key_here

4. Run scripts step by step:
   > python 1.fetch.py  
   > python 2.tokenization.py  
   > python 3.setup_qdrant.py  
   > python 4.embedding_insert.py  

5. Launch Streamlit app:
   > streamlit run app.py

-------------------------------------
Query Example
-------------

Question:
What are the treatment options for pancreatic cancer in elderly patients?

- Fetches top 5 most relevant papers
- Shows titles, links, and scores
- Gemini provides a final summarized answer

-------------------------------------
Notes
-----

- Everything runs locally, no Docker or external DBs needed
- You can extend this to other domains by changing the dataset
- Gemini output is based only on retrieved content

-------------------------------------
License
-------

This project is for research and educational purposes only.
Not intended for clinical use.
