"""Search PubMed and return a pandas DataFrame of recent articles."""
from __future__ import annotations

import os
from datetime import date, timedelta
from typing import List, Dict

import pandas as pd
from Bio import Entrez
from xml.etree import ElementTree as ET

# Use the email you set in your GitHub secret EMAIL_USER
Entrez.email = os.getenv("EMAIL_USER", "you@example.com")


def _xml_to_dict(article_xml) -> Dict:
    """Extract fields from a <PubmedArticle> element into a dict."""
    medline = article_xml.find("MedlineCitation")
    article = medline.find("Article")

    journal = article.find("Journal").find("Title").text
    title = article.find("ArticleTitle").text

    # Concatenate all <AbstractText> paragraphs
    abstract_tag = article.find("Abstract")
    abstract = (
        " ".join(p.text for p in abstract_tag.findall("AbstractText"))
        if abstract_tag is not None else ""
    )

    pmid = medline.find("PMID").text

    # Find DOI if present
    doi = None
    for id_tag in article.findall("ELocationID"):
        if id_tag.attrib.get("EIdType") == "doi":
            doi = id_tag.text
            break

    # Build author list (up to 6 + “et al.”)
    authors = []
    for author in article.findall("AuthorList/Author"):
        last = author.findtext("LastName", "")
        first = author.findtext("ForeName", "")
        if last:
            authors.append(f"{last} {first[0]}." if first else last)
    authors_str = ", ".join(authors[:6]) + (" et al." if len(authors) > 6 else "")

    return {
        "pmid": pmid,
        "title": title,
        "authors": authors_str,
        "journal": journal,
        "doi": doi,
        "abstract": abstract,
    }


def fetch_last_week(
    keywords: List[str],
    window_days: int = 7,
    retmax: int = 100,
    api_key: str | None = None
) -> pd.DataFrame:
    """
    Query PubMed for articles matching `keywords` published in the last `window_days` days.
    Returns a DataFrame with columns: pmid, title, authors, journal, doi, abstract.
    """
    term = " OR ".join(keywords)
    today = date.today()
    mindate = today - timedelta(days=window_days)

    # ESearch to get PMIDs
    search_handle = Entrez.esearch(
        db="pubmed",
        term=term,
        mindate=mindate.strftime("%Y/%m/%d"),
        maxdate=today.strftime("%Y/%m/%d"),
        retmax=retmax,
        api_key=api_key,
    )
    id_list = Entrez.read(search_handle)["IdList"]

    # No results → empty DataFrame
    if not id_list:
        return pd.DataFrame()

    # EFetch to get full XML for those PMIDs
    fetch_handle = Entrez.efetch(
        db="pubmed",
        id=",".join(id_list),
        rettype="xml",
        api_key=api_key,
    )
    tree = ET.parse(fetch_handle)
    root = tree.getroot()

    # Convert each <PubmedArticle> into a dict row
    rows = [_xml_to_dict(article) for article in root.findall("PubmedArticle")]

    return pd.DataFrame(rows)
