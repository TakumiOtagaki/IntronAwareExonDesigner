# Evaluator design notes

- Step 1 focuses on translating amino-acid sequences into codons with Biopython's NCBI table to keep the decision space purely synonymous.
- ViennaRNA is imported lazily so that the module can still be inspected even if the bindings are unavailable; evaluation errors are surfaced with a clear runtime message.
- Boundary pair scoring aggregates probabilities for the intron termini that the context reports, while caching avoids repeated ViennaRNA runs for identical sequences.
