"""Send email with SMTP."""
import os
import smtplib
from email.message import EmailMessage
from pathlib import Path

SMTP_HOST = "smtp.mail.me.com"
SMTP_PORT = 587


def send_mail(subject: str, body: str, attachment: Path | None = None,
              to_addr: str | None = None) -> None:
    user = os.getenv("EMAIL_USER")
    pwd = os.getenv("EMAIL_PASS")
    to_addr = to_addr or user

    msg = EmailMessage()
    msg["From"] = user
    msg["To"] = to_addr
    msg["Subject"] = subject
    msg.set_content(body)

    if attachment and attachment.exists():
        msg.add_attachment(attachment.read_bytes(), maintype="text", subtype="markdown",
                           filename=attachment.name)

    with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as smtp:
        smtp.starttls()
        smtp.login(user, pwd)
        smtp.send_message(msg)