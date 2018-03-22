import dawg

# akita = dawg.DawgWalker(b"../examples/basic-dna.dawg", b"fasta:-", 10, 121212)
akita = dawg.DawgWalker()
# akita.bark(b"bark bark bark")
akita.run(b"../examples/basic-dna.dawg", b"fasta:-", 10, 212121)
