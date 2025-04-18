# PubMed Weekly Scraper

Automated GitHub Action that every Friday:

1. Searches PubMed for the past 7 days using your custom keywords.
2. Retrieves titles, authors, journals, DOIs, and abstracts.
3. Summarises each abstract with OpenAI into a 3‑sentence "Why this matters" blurb.
4. Builds a Markdown report **and** e‑mails it to you.

## ⚡ Quick start (GUI‑only)

1. Click **Use this template → Create a new repo** *(or clone if you already downloaded the zip)*.
2. Open the repo in **GitHub Desktop** and **Clone** to your Mac.
3. In VS Code, open the folder and paste the files below into the correct paths.
4. **Commit → Push**.
5. In **Repo → Settings → Secrets → Actions** add:
   * `OPENAI_API_KEY`  →  `sk-…`
   * `EMAIL_USER`      →  your sending Gmail
   * `EMAIL_PASS`      →  App Password
   * `NCBI_KEY`        →  *(optional)*
6. **Actions** tab → select **weekly** → **Run workflow**.
7. Check your inbox 🎉  (or see `output/` folder created in the repo).

That’s it—every Friday 13:00 ET a fresh report appears.