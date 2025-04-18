"""Render Markdown report with Jinja2."""
import datetime as _dt
from pathlib import Path
from typing import Iterable

import jinja2
import pandas as pd

_TEMPLATE = """## PubMed Ortho & Spine Scan – Week of {{ today }}

{% for row in rows %}
### {{ row.title }}
*{{ row.authors }} – *{{ row.journal }}*

{{ row.summary }}
{% if row.doi %}
[PubMed](https://pubmed.ncbi.nlm.nih.gov/{{ row.pmid }}/) | [DOI](https://doi.org/{{ row.doi }})
{% else %}
[PubMed](https://pubmed.ncbi.nlm.nih.gov/{{ row.pmid }}/)
{% endif %}

---
{% endfor %}
"""


def write_markdown(df: pd.DataFrame, out_dir: Path = Path("output")) -> Path:
    out_dir.mkdir(exist_ok=True)
    env = jinja2.Environment(trim_blocks=True, lstrip_blocks=True)
    template = env.from_string(_TEMPLATE)
    md = template.render(today=_dt.date.today().strftime("%d %b %Y"), rows=df.to_dict("records"))
    path = out_dir / f"pubmed_{_dt.date.today().isoformat()}.md"
    path.write_text(md, encoding="utf-8")
    return path