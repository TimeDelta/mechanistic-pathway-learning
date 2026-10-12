"""The limits of the shell fallback for public pages (experiments/read_public_page.py), checked without the network."""
from __future__ import annotations

import importlib.util
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parents[1] / "experiments" / "read_public_page.py"
module_specification = importlib.util.spec_from_file_location("read_public_page", MODULE_PATH)
read_public_page = importlib.util.module_from_spec(module_specification)
module_specification.loader.exec_module(read_public_page)

ROBOTS_FILE_THAT_TURNS_ASSISTANTS_AWAY = """
User-agent: *
Disallow: /private/

User-agent: ClaudeBot
Disallow: /
"""

ROBOTS_FILE_WITH_ONE_CLOSED_FOLDER = """
User-agent: *
Disallow: /private/
"""


def test_an_address_with_an_email_address_is_not_requested():
    assert read_public_page.reason_address_is_not_requested("https://api.example.org/v2/10.1/x?email=someone%40example.org") is not None
    assert read_public_page.reason_address_is_not_requested("https://api.example.org/v2/10.1/x?email=someone@example.org") is not None


def test_only_public_https_addresses_without_credentials_are_requested():
    assert read_public_page.reason_address_is_not_requested("http://example.org/page") is not None
    assert read_public_page.reason_address_is_not_requested("https://localhost/page") is not None
    assert read_public_page.reason_address_is_not_requested("https://user:secret@example.org/page") is not None
    assert read_public_page.reason_address_is_not_requested("https://example.org/page?id=7") is None


def test_a_rule_against_the_assistant_crawler_closes_the_whole_site():
    assert not read_public_page.robots_file_allows(ROBOTS_FILE_THAT_TURNS_ASSISTANTS_AWAY, "https://example.org/open/page")


def test_a_rule_for_every_reader_closes_only_what_it_names():
    assert read_public_page.robots_file_allows(ROBOTS_FILE_WITH_ONE_CLOSED_FOLDER, "https://example.org/open/page")
    assert not read_public_page.robots_file_allows(ROBOTS_FILE_WITH_ONE_CLOSED_FOLDER, "https://example.org/private/page")


def test_an_empty_robots_file_allows_every_address():
    assert read_public_page.robots_file_allows("", "https://example.org/any/page")


def test_a_web_page_served_in_place_of_the_robots_file_is_not_a_robots_file():
    assert not read_public_page.looks_like_robots_file("<!DOCTYPE html><html><head><title>Verify</title></head></html>")
    assert read_public_page.looks_like_robots_file("User-agent: *\nDisallow:\n")


def test_a_human_check_is_recognised_and_an_ordinary_page_is_not():
    assert read_public_page.looks_like_challenge(403, "Forbidden")
    assert read_public_page.looks_like_challenge(200, "<html><body>To continue, we need to verify you are a human.</body></html>")
    # a proof-of-work gate answers with status 200 and a short page of its own
    assert read_public_page.looks_like_challenge(200, "<html><head><title>Making sure you're not a bot!</title></head><body><script src='/.within.website/x/cmd/anubis/static/js/main.mjs'></script></body></html>")
    assert read_public_page.looks_like_challenge(200, "<html><body>Checking your browser before accessing the site.</body></html>")
    long_article_that_mentions_a_captcha = "<html><body>" + "A study of captcha design. " * 2000 + "</body></html>"
    assert not read_public_page.looks_like_challenge(200, long_article_that_mentions_a_captcha)
    assert not read_public_page.looks_like_challenge(200, "<html><body>Release 24.0, 16 March 2026.</body></html>")


def test_a_refused_address_makes_no_request(tmp_path):
    reader = read_public_page.PublicPageReader(tmp_path)
    record = reader.read("http://example.org/page")
    assert record["read"] is False and "https" in record["reason"]
    assert not any(tmp_path.iterdir())


def test_text_of_page_drops_scripts_and_tags():
    page_source = "<html><head><style>p {color: red}</style><script>var a = 1;</script></head><body><p>Dose &amp; effect</p></body></html>"
    assert read_public_page.text_of_page(page_source) == "Dose & effect"


def test_a_proof_of_work_gate_whose_challenge_is_only_in_its_script_is_a_refusal():
    """PubMed Central's scanned-PDF gate answers 200 with a short interstitial whose only challenge markers are in its
    script, so the prose markers missed it and a research pass of 12 October 2026 saved it as if it were the document."""
    interstitial = ('<!DOCTYPE html><html><head><title>Preparing to download</title></head><body>'
                    '<p>Preparing to download ...</p><script>var POW_CHALLENGE="ab12"; var POW_DIFFICULTY=4;'
                    'window.ncbi.pmc.pow.init(POW_CHALLENGE, POW_DIFFICULTY);</script></body></html>')
    assert read_public_page.looks_like_challenge(200, interstitial)
    assert not read_public_page.looks_like_challenge(200, "<!DOCTYPE html><html><body><p>" + "the article itself. " * 2000 + "</p></body></html>")
