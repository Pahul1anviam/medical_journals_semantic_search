import json
from nltk.tokenize import sent_tokenize
import os
from nltk.tokenize.punkt import PunktSentenceTokenizer
import nltk

nltk.download("punkt")

def chunk_text(text, max_tokens=250):
    tokenizer = PunktSentenceTokenizer()
    sentences = tokenizer.tokenize(text)
    chunks, current_chunk, current_len = [], "", 0

    for sentence in sentences:
        sentence_len = len(sentence.split())
        if current_len + sentence_len > max_tokens:
            chunks.append(current_chunk.strip())
            current_chunk = sentence
            current_len = sentence_len
        else:
            current_chunk += " " + sentence
            current_len += sentence_len
    if current_chunk.strip():
        chunks.append(current_chunk.strip())
    return chunks

def preprocess_abstracts(input_file, output_file):
    with open(input_file, "r", encoding="utf-8") as f_in, open(output_file, "w", encoding="utf-8") as f_out:
        for line in f_in:
            data = json.loads(line)
            abstract = data.get("abstract", "")
            if len(abstract.strip()) > 100:
                chunks = chunk_text(abstract)
                for chunk in chunks:
                    item = {
                        "chunk": chunk,
                        "pmid": data.get("pmid", ""),
                        "title": data.get("title", ""),
                        "authors": data.get("authors", []),
                        "url": data.get("url", ""),
                        "journal": data.get("journal", "")
                    }
                    f_out.write(json.dumps(item) + "\n")
    print("✅ Tokenized and saved with metadata")

if __name__ == "__main__":
    preprocess_abstracts("pubmed_abstracts_new.jsonl", "pubmed_chunks_new.jsonl")
