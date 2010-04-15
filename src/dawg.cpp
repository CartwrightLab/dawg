/*  Dawg - DNA Assembly with Gaps - Simulating Sequence Evolution
    Copyright (c) 2004-2009  Reed A. Cartwright, PhD

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <boost/preprocessor.hpp>
#include <boost/config.hpp>

#include <exception>

#include "dawg.h"
#include "rand.h"
#include "var.h"
#include "dawg_app.h"

#include <dawg/ma.h>
#include <dawg/pile.h>
#include <dawg/utils/foreach.h>
#include <dawg/matic.h>
#include <dawg/global.h>

using namespace std;
using namespace boost;
using namespace dawg;

#define VERSION_MSG PACKAGE_STRING "\n" \
	"    Copyright (C) 2005-2009  Reed A. Cartwright, PhD <reed@scit.us>\n" \
	"    Report bugs to " PACKAGE_BUGREPORT

// Help Information
char g_csDawgTxt[] =
#include "dawgtxt.h"
;

bool Parse(const char* cs);
bool Execute();

bool g_bReportErrors = true;
bool g_bReportWarnings = true;

const char *g_csOutput = NULL;

int main(int argc, char *argv[])
{
	int ret = EXIT_FAILURE;
	try {
		dawg_app app(argc, argv);
		ret = app.run();
	} catch(std::exception &e) {
		CERROR(e.what());
	}
	return ret;
}

dawg_app::dawg_app(int argc, char* argv[]) : desc("Allowed Options") {
	try {
		desc.add_options()
			#define XM(lname, sname, desc, type, def) ( \
				_S(lname) _IFD(sname, "," BOOST_PP_STRINGIZE sname), \
				po::value< type >(&arg._V(lname))->default_value(def), \
				desc )				
			#include "dawgarg.xmh"
			#undef XM
			;
		po::variables_map vm;
		po::positional_options_description pdesc;
		pdesc.add("input", -1);
		po::store(po::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vm);
		po::notify(vm);
		if(!arg.arg_file.empty()) {
			if(arg.arg_file == "-") {
				po::store(po::parse_config_file(cin, desc), vm);	
			} else {
				std::ifstream ifs(arg.arg_file.c_str());
				if(!ifs.is_open()) {
					string sse = "unable to open argument file ";
					sse += arg.arg_file;
					throw std::runtime_error(sse);
				}
				po::store(po::parse_config_file(ifs, desc), vm);
			}
			po::notify(vm);
		}
	} catch (std::exception &e) {
		CERROR(e.what());
		throw std::runtime_error("unable to process command line");
	}
}

int dawg_app::run() {
	if(arg.version)	{
		cerr << endl << VERSION_MSG << endl << endl;
		return EXIT_SUCCESS;
	}
	if(arg.help) {
		cerr << endl << VERSION_MSG << endl << endl;
		cerr << desc << endl;
		return EXIT_SUCCESS;
	}
	//if(arg.quiet)
	//	cerr.clear(ios::failbit);
	pile input;
		
	bool ret = true;
	foreach(string &ss, arg.input) {
		ret &= input.parse_file(ss.c_str());
	}
	if(!ret)
		return EXIT_FAILURE;	
	// process aliases
	input.read_aliases();

	global_options glopts;
	glopts.read_section(input.data.front());

	vector<dawg::ma> configs;
	if(!dawg::ma::from_pile(input, configs)) {
		DAWG_ERROR("bad configuration");
		return EXIT_FAILURE;
	}
		
	// Create the object that will do all the simulation
	// work for us.  Configure its sections.
	dawg::matic kimura;
	if(!kimura.configure(configs.begin(), configs.end())) {
		DAWG_ERROR("bad configuration");
		return EXIT_FAILURE;
	}
	// if a seed was specified, use it
	if(!glopts.sim_seed.empty()) {
		kimura.seed(glopts.sim_seed.begin(), glopts.sim_seed.end());
	}
	// create sets of aligned sequences;
	dawg::alignment aln;
	for(unsigned int i=1;i<=glopts.sim_reps;++i) {
		kimura.walk(aln);
		cout << aln << endl;
	}
		
	return EXIT_SUCCESS;
}

// Fast and dirty random seed algorithm
// based off of NR's randqd1
inline unsigned int rand_seed()
{
	static unsigned int u = (unsigned int)(time(NULL)+3*getpid());
	return (u = u*1664525u + 1013904223u);
}

inline const char * SS2CS(const string& ss)
{
	return ss.empty() ? NULL : ss.c_str();
}

// Execute the Dawg process based on the current configuration
bool Execute()
{
	// Variables
	unsigned int uReps  =  1;
	vector<unsigned int> vuSeqLen;
	vector<string> vssSeqs;
	string::size_type uTotalSeqLen = 0, uTotalRateLen = 0;

	double dGamma, dIota;

	double dNucFreq[4] = {0.25,0.25,0.25,0.25};
	double dRevParams[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	vector<double> vdParams;
	vector< vector<double> > vvdRates;
	string ssModel = "JC", ssGapModel[2] = {"US", "US"};

	double dLambda[2] = {0.0, 0.0};
	vector<double> vdGapModel[2];
	vdGapModel[0].push_back(1.0);
	vdGapModel[1].push_back(1.0);
	unsigned int uKeepFlank = 0;

	vector<NewickNode*> vtTrees;
	double	dTreeScale = 1.0;
	vector<unsigned int>   vuSeed;

	string ssFile = "-";
	unsigned int uFmt = FormatClustal;
	string ssFormat = "";
	string ssOutBlockHead = "";
	string ssOutBlockBefore = "";
	string ssOutBlockAfter = "";
	string ssOutBlockTail = "";
	bool bOutSubst = true;

	bool bGapSingle = false, bGapPlus = false, bLowerCase = false, bTranslate = false;
	bool bKeepEmpty = false;
	DawgVar::Vec::size_type nRes;

	// Read variables from configuration


	// Read Out.Block.* from file if neccessary
   	if(!ssOutBlockHead.empty())
	{
		ifstream iFile(ssOutBlockHead.c_str());
		if(iFile.is_open())
			getline(iFile, ssOutBlockHead, '\0');
	}
   	if(!ssOutBlockBefore.empty())
	{
		ifstream iFile(ssOutBlockBefore.c_str());
		if(iFile.is_open())
			getline(iFile, ssOutBlockBefore, '\0');
	}
   	if(!ssOutBlockAfter.empty())
	{
		ifstream iFile(ssOutBlockAfter.c_str());
		if(iFile.is_open())
			getline(iFile, ssOutBlockAfter, '\0');
	}
   	if(!ssOutBlockTail.empty())
	{
		ifstream iFile(ssOutBlockTail.c_str());
		if(iFile.is_open())
			getline(iFile, ssOutBlockTail, '\0');
	}

	// Output Format
	if(ssFormat == "")
	{
		string::size_type pos = ssFile.find_first_of(':');
		if( pos != string::npos)
			ssFormat = ssFile.substr(0, pos);
		else
		{
			pos = ssFile.find_last_of('.');
			if(pos != string::npos)
				ssFormat = ssFile.substr(pos+1);
			else
				ssFormat = "Clustal";
		}
	}
	// Override output
	if(g_csOutput != NULL)
	{
		ssFile = g_csOutput;
		string::size_type pos = ssFile.find_first_of(':');
		if( pos != string::npos)
			ssFormat = ssFile.substr(0, pos);
		else
		{
			pos = ssFile.find_last_of('.');
			if(pos != string::npos)
				ssFormat = ssFile.substr(pos+1);
		}
	}
	if(ssFormat == "Fasta" || ssFormat == "fas")
		uFmt = FormatFasta;
	else if(ssFormat == "Nexus" || ssFormat == "nex")
		uFmt = FormatNexus;
	else if(ssFormat == "Phylip" || ssFormat == "phy")
		uFmt = FormatPhylip;
	else if(ssFormat == "Clustal" || ssFormat == "aln"
		|| ssFormat == "poo" || ssFormat == "txt" || ssFormat == "out")
		uFmt = FormatClustal;
	else
		return DawgError("Unknown file format, \"%s\".", ssFormat.c_str());

	// setup output location
	ostream* pOut;
	ofstream ofOut;
	string::size_type pos = ssFile.find_first_of(':');
	if( pos != string::npos)
		ssFile.erase(0, pos+1);

	if(ssFile == "-" || ssFile.empty())
		pOut = &cout;
	else
	{
		ofOut.open(ssFile.c_str());
		if(!ofOut.is_open())
			return DawgError("Unable to open \"%s\" for output.", ssFile.c_str());
		pOut = &ofOut;
	}
	DawgIniOutput(*pOut);

	// Evolve Many Sequences
	return true;
}

// Error Reporting Routine
bool DawgError(const char* csErr, ...)
{
	if(!g_bReportErrors)
		return false;
	fprintf(stderr, "Error: ");
	va_list args;
	va_start(args, csErr);
#if defined(HAVE_VFPRINTF)
	vfprintf(stderr, csErr, args);
#elif defined(HAVE__DOPRNT)
	_doprnt(csErr, args, stderr);
#else
	fprintf(stderr, "%s", csErr);
#endif
	va_end(args);
	fprintf(stderr, "\n");
	return false;
}

bool DawgWarn(const char* csErr, ...)
{
	if(!g_bReportWarnings)
		return false;
	fprintf(stderr, "Warning: ");
	va_list args;
	va_start(args, csErr);
#if defined(HAVE_VFPRINTF)
	vfprintf(stderr, csErr, args);
#elif defined(HAVE__DOPRNT)
	_doprnt(csErr, args, stderr);
#else
	fprintf(stderr, "%s", csErr);
#endif
	va_end(args);
	fprintf(stderr, "\n");
	return false;
}
