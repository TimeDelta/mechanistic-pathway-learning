# Working rules for sessions in this repository

These are the owner's rulings on how sources are read. Each is dated and quoted. They hold for every session and
for every research pass a session starts, so a session passes them on to the passes it starts.

## Reading sources (11 October 2026)

- **An unanswered permission request is not a finding.** "Unacceptable to withdraw due to lack of answered
  permission request." The fetch tool withdraws a request after five minutes without an answer. Keep a list of the
  requests that ended that way, repeat each one as soon as the owner is back, tell the owner at once which pages
  were affected, and never write "not retrieved" for that reason.
- **Why the request appears, and how to avoid needing it.** In the Claude app a fetch "is limited to search results
  and URLs you've shared" (support.claude.com, "Get started with Claude Cowork", read 11 October 2026). An address
  that came from neither needs the owner's approval. Neither the `WebFetch` allow rule in `.claude/settings.json` nor the
  app's Auto mode removed the prompt in a test on 11 October 2026. So: search first, so that the address is a search result; use the
  connectors and the providers' programmatic services where they exist; and before work that will run unattended,
  send the owner the list of sites the work will need.
- **A licence is verified.** "No. You have to verify the license." When the provider's licence page cannot be
  read, read the terms where the provider also states them: the header of the data file, the provider's source
  repository, the data availability statement of the provider's paper. Quote the sentence and say where it was
  read. When every such route is closed, name the licence as unverified, list the routes tried and give the address
  for the owner to open.
- **A refusal by the site is respected.** A page behind a robots rule or a check that the visitor is human is not
  read through a cache, a mirror or a solved challenge. Record the refusal and use the provider's other
  publications.

The table of licences is `docs/data_sources.md`. The quotes and where each was read are in
`docs/label_source_review.md`, section 11.
