# An Analysis of Dawg as a Python Module
 - The goal of these changes is to use DAWG effectively in a high-level language like Python
 - Use DAWG in Python should make it easier to write tests and integrate with other Python libraries
 - Boost.Python, Cython, and SWIG are the three main wrapping tools analyzed
 - A BASH script was used to do the time trials **analyze_wrappers.sh**
    - Each wrapper was tested with the **basic-dna.dawg** trick file
    - Each wrapper used 10 reps and seed **212121**
    - Each wrapper program was called continuously 10000 times, output redirected to **/dev/null**
    - Linux **time** command was used to record the time (real time)
 - The sizes of their libraries were analyzed with **ls -lh**
 - A simple API was created to try to show how to use DAWG quickly from Python
    - This API was more or less consistent between all three wrappers
    - It consisted of a **run** command which took a DAWG trick file, output, rep, and seed
    - There are some discrepencies between the APIs and how they are called:
        - For instance, Cython requires that a string be converted into bytes before getting called

```c++
void dawg::Dawg::run(
    const std::string &inFile,
    const std::string &outFile,
    const unsigned int reps,
    const unsigned int seed)
{
    using std::string;

	dawg::trick input;

    bool ret &= dawg::trick::parse_file(input, inFile.c_str());

	if(!ret)
		std::cerr << "Failure to parse DAWG file\n";

	input.read_aliases();

	dawg::global_options glopts;
	glopts.read_section(input.data.front());

	dawg::output write_aln;

	if (!write_aln.open(/*glopts.output_file.c_str()*/ outFile.c_str(),
		reps - 1,
		false,
		false,
		false))
	{
		DAWG_ERROR("bad configuration");
		return;
	}
	write_aln.set_blocks(glopts.output_block_head.c_str(),
		glopts.output_block_between.c_str(),
		glopts.output_block_tail.c_str(),
		glopts.output_block_before.c_str(),
		glopts.output_block_after.c_str()
	);

	std::vector<dawg::ma> configs;
	if (!dawg::ma::from_trick(input, configs)) { // throwing error
		DAWG_ERROR("bad configuration");
		return;
	}

	dawg::matic kimura;
	if(seed != 0) {
		kimura.seed(this->seed);
	}

	if (!kimura.configure(configs.begin(), configs.end())) {
		DAWG_ERROR("bad configuration");
		return;
	}

	// create sets of aligned sequences;
    std::vector<dawg::alignment> alignments;
	dawg::alignment aln;
	kimura.pre_walk(aln);
	for (unsigned int i = 0; i< reps; ++i) {
		kimura.walk(aln);
		alignments.insert(alignments.end(), aln);
		write_aln(aln); // this would print the aln data out to std::cout or a file
    }
} // run
```

 - A builder-style call chain can be created as well (not implemented)
```python
alignments = PyDawg.trick("").walk("").fetch("")
for i in range(alignments):
    print("alignment[{}]: {}".format(i, alignments[i]))
```

    - The RNG within DAWG is easy to use:
```c++
unsigned int
dawg::Dawg::rand(unsigned int a, unsigned int b) {
    static dawg::mutt rng;
    //rng.seed(11);
    return rng.rand_uint() % b + a;
}
```
```python
rng = PyDawg.PyDawg()
rng.rand(0, 100)
```

## Boost.Python
 - Time trial: 6m45.847s
 - **dawg.so** size: 6.7M

## Cython
 - Time trial: 3m51.320s
 - **PyDawg.cpython-36m-x86_64-linux-gnu.so** size: 20M

## SWIG
 - Time trial: 7m17.167
 - **_dawg.so** size: 6.6M


 # Additional Notes
  - Working on an example that uses Python methods to parse DAWG file
  and map directly to the DAWG trick/matic classes
  - Could try integrating with SISRS?



### General Python Notes
 - list comprehension syntax: [ expression for item in list if conditional ]

 - generator expressions are functions that yield their iterators ...



```python
# The original example I got working
# Must use bytes on the string from Cython -> CPP (Cython book)
# donovan = pd.PyDawg(b"../examples/basic-dna.dawg", b"fasta:-", 10, 212121)
# donovan.run()

# PyDawg constructor can take in Tricks, Walks, and Bones
# akita = pd.pd(
# DawgBuilder(
#     Trick(name="BasicDnaExample",
#         Segments(name="__default__", inheritsFrom="LUCA",
#             (tree_tree="((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);",
#             root_length=1000, subst_model='hky',
#             subst_freqs=[0.2, 0.3, 0.3, 0.2],
#             subst_params=[0.2, 1.0], sim_reps=10))),
#     Walk(), Bone("fasta:-")))
```
