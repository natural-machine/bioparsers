"""Network retrieval of reference-database entries by identifier.

The complement to ``bioparsers.parsers``: where a parser reads a local
distribution file top to bottom, a fetcher pulls the specific entries named by
a list of identifiers. Both hand back the same ``Record`` types, so anything
downstream — builders, recipes, ``dump_jsonl`` — cannot tell which was used.

Reach for a fetcher when the wanted subset is small relative to the
distribution file. Retrieving a few thousand TrEMBL entries over the REST API
takes seconds; finding them by scanning the 160 GB local
``uniprot_trembl.dat.gz`` takes hours. The trade is currency and reachability:
a fetcher returns whatever UniProt currently serves (so accessions merged or
deleted since a local snapshot will not match it, and the annotation may have
drifted), and it needs the network.

One module per source. Currently ``uniprot_rest`` (UniProtKB).
"""
