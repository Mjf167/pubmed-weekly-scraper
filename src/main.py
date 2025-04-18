"""Entry‑point orchestrator."""
import os
import yaml

import pandas as pd

from pathlib import Path

from query_pubmed import fetch_last_week
from summarize import summarise
from render import write_markdown
from send_email import send_mail

CONFIG_PATH = Path(__file__).parent.parent / "config.yml"


def load_cfg():
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main() -> None:
    cfg = load_cfg()
    df = fetch_last_week(cfg["keywords"], cfg["window_days"], cfg["max_papers"], os.getenv("NCBI_KEY"))

    if df.empty:
        print("No new papers this week.")
        return

    # Add summaries
    df["summary"] = df["abstract"].apply(lambda abs_: summarise(abs_, cfg["summary_model"]))

    md_path = write_markdown(df)
    send_mail(subject="Weekly Ortho PubMed Scan", body="Report attached.", attachment=md_path,
              to_addr=cfg["email_to"])


if __name__ == "__main__":
    main()