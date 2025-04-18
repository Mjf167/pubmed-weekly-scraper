"""OpenAI summarisation helper."""
import os
from openai import OpenAI

SYSTEM_PROMPT = (
    "You are an expert orthopaedic surgery educator. "
    "Write a concise, three‑sentence summary for an MS4 interested in spine surgery. "
    "Finish with one sentence: 'Why it matters: …'."
)


def summarise(text: str, model: str = "gpt-4o") -> str:
    client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))
    resp = client.chat.completions.create(
        model=model,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": text},
        ],
        max_tokens=120,
    )
    return resp.choices[0].message.content.strip()