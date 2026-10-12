# Notice of use to CTD: a draft for the owner to send

CTD's data files state: "You must notify CTD, and describe your use of our data" and "For quality control purposes,
you must provide CTD with periodic access to your publication of our data." Only the repository's owner can send
that notice. The form is at <http://ctdbase.org/help/contact.go> (the site declines automated readers, so the form's
fields were not seen; the text below is written to be pasted into a free-text field).

When it is sent, record the date in the CTD row of `docs/data_sources.md` and in the decision log.

Before sending, two things are worth the owner's own reading, because neither page could be read from here:

- the legal notice, <http://ctdbase.org/about/legal.jsp>. The data file says its terms "include" the four it lists, so
  there may be others. The point to look for is whether a table derived from CTD may be redistributed in a public
  repository. The draft asks CTD that question directly;
- the citation guidelines, <http://ctdbase.org/about/publications/#citing>, in case they ask for wording other than
  the citation of the 2025 paper that the repository carries.

## Draft

Subject: Notice of use of CTD chemical-disease data in a public research repository

Dear CTD team,

I am writing to notify you of my use of CTD data, as your terms ask.

Who: [name, affiliation or "independent researcher", contact address]

Project: "mechanistic-pathway-learning", a non-commercial research study that predicts psychiatric and neurological
symptoms from gene and drug perturbations with a graph model of physiology. The repository is public:
https://github.com/TimeDelta/mechanistic-pathway-learning

Data used: the file CTD_chemicals_diseases.tsv.gz (report created Tue Sep 29 13:17:30 EDT 2026, downloaded on
3 October 2026). From it I keep the rows with direct evidence ("marker/mechanism" or "therapeutic") whose disease is one
of the MeSH descriptors of the study's 24 symptoms and whose chemical matches a drug of the study. That gives 5,191
rows, one per cited paper, from 3,356 papers and 457 drugs. Rows inferred through a gene are not used.

How they are used: as literature evidence of low grade beside drug-label and genetic evidence. They inform a
reliability model as a soft prior. They are never used as labels for evaluation.

Publication of the data: the derived table is in the repository at
data/releases/v0.4/literature/ctd_relations_filtered.csv, with the CTD chemical and disease identifiers and the PubMed
identifier of each row. Beside it, the file CTD_TERMS.md quotes your terms, cites Davis et al. 2025
(doi:10.1093/nar/gkae883), states that the data are for non-commercial use and gives the form of the links to your
chemical and disease pages. The repository's README carries the citation as well. The repository is public, which I
hope meets your request for access; I am glad to give access in another form.

Two questions:

1. Is it acceptable to you that this derived table is redistributed in a public repository, with your terms and the
   citation beside it? If not, I will remove it and have the build script download the file from CTD instead.
2. Is there wording for the citation or for the links that you prefer to what the repository carries now?

Thank you for maintaining CTD and for making it available to researchers.

[name]
