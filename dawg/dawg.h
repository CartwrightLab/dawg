#ifndef DAWG_DAWG_H
#define DAWG_DAWG_H

#include <vector>
#include <string>
#include <fstream>

#include "node.h"

inline const char* DawgVersion()
{
	return "1.1.0";
}

// Error Reporting
bool DawgError(const char* csErr, ...);  //always returns false

// Functions to set up evolution model
bool EvoRevParams(double pFreqs[], double pSubs[]);
bool EvoRateParams(double dGamma, double dIota);
bool EvoScaleTrees(double dScale);
bool EvoIndelUser(double dInsL, double pIns[], int nIns,double dDelL, double pDel[], int nDel);
bool EvoIndelNegBn(double dInsL, unsigned long uInsR, double dInsQ, double dDelL, unsigned long uDelR, double dDelQ);

std::string EvoDescription();

// Sequence Output
class MapSsToVar;
class DawgVar;
class Node;
typedef Node Tree;

enum FileFormat { FASTA, NEXUS, PHYLIP };
bool SetFormat(FileFormat fmt, int nNum, const char* csBlock, bool bGapSingle);
bool DawgOpen(const char* csFile, std::ofstream& rFile);
void DawgIniOutput(std::ostream& os);
bool SaveSequences(std::ostream &rFile, const Tree *arTrees[], unsigned int uSize);

class DawgVar
{
public:
	typedef std::vector<DawgVar*> Vec;

	enum Type { tyNone, tyBool, tyNumber, tyString, tyVector, tyTree};

	DawgVar() : m_tyType(tyNone) { }
	virtual ~DawgVar();
	
	// Type Routines
	Type GetType() const { return m_tyType; }
	bool IsType(Type ty) const { return m_tyType == ty; }
	void Unset();

	// Variable Retreival / Setting
	static MapSsToVar& GetMap();
	static DawgVar* GetVar(const std::string &ssKey);
	static void		SetVar(const std::string &ssKey, DawgVar* pVar);
	static void		ClearMap();
	
	// General Routines
	unsigned int Size();

	// Number Routines
	explicit DawgVar(double dVar) : m_tyType(tyNumber), m_dData(dVar) { }
	double		GetNumber() const { return m_dData; }
	bool		Get(double &rdVar ) const;
	void		Set(double dVar );
	bool		Get(int &rnVar ) const;
	void		Set(int nVar );
	
	// Bool Routines
	explicit DawgVar(bool bVar) : m_tyType(tyBool), m_bData(bVar) { }
	bool		GetBool() const { return m_bData; }
	bool		Get(bool &rbVar ) const;
	void		Set(bool bVar );

    // String Routines
	explicit DawgVar(const std::string &rssVar) : m_tyType(tyString)
		{ m_pssData = new std::string(rssVar); }
	const std::string&	GetString() const { return *m_pssData; }
	bool				Get(std::string &rssVar ) const;
	void				Set(const std::string &rssVar );
	
	// Vector Routines
	explicit DawgVar(const Vec *p) : m_tyType(tyVector)
		{ m_pvData = p; }
	const Vec*	GetVector() const { return m_pvData; }
	bool		Get(const Vec *&rVec) const;
	void		Set(const Vec* rVec);

	DawgVar& GetAt(Vec::size_type uIndex)
		{ return IsType(tyVector) ? *(*m_pvData)[uIndex] : *this; }
	const DawgVar& GetAt(Vec::size_type uIndex) const
		{ return IsType(tyVector) ? *(*m_pvData)[uIndex] : *this; }
	DawgVar& operator[](Vec::size_type uIndex)
		{ return GetAt(uIndex); }
	const DawgVar& operator[](Vec::size_type uIndex) const
		{ return GetAt(uIndex); }
	
	// Tree Routines
	explicit DawgVar(Tree* p) : m_tyType(tyTree), m_ptrData(p) { }
	Tree*	GetTree() const { return m_ptrData; }
	bool		Get(Tree *&rTree) const;
	void		Set(Tree *pTree);

protected:
	Type m_tyType;
	union
	{
		double m_dData;
		bool   m_bData;
		const Vec*   m_pvData;
		const std::string* m_pssData;
		Tree*  m_ptrData;
	};
private:
	DawgVar(const DawgVar& var);
	DawgVar& operator = (const DawgVar& var);

public:
	// Templates
	template<class _T>
	static bool Get(const std::string &ssKey, _T &r)
	{
		DawgVar *pVar = GetVar(ssKey);
		return ( pVar != NULL && pVar->Get(r) );
	}
	template<class _T>
	Vec::size_type GetArray(_T ar[], Vec::size_type uSize)
	{
		Vec::size_type uMax = min(uSize, Size());
		Vec::size_type u = 0;
		for(; u<uMax; u++)
		{
			if(!GetAt(u).Get(ar[u]))
				return u;
		}
		return u;		
	}

	template<class _T>
	bool GetVector(std::vector<_T> &rVec)
	{
		_T tTemp;
		rVec.clear();
		for(Vec::size_type u = 0; u<Size(); u++)
		{
			if(!GetAt(u).Get(tTemp))
				return false;
			rVec.push_back(tTemp);
		}
		return true;
	}
	template< class _T >
	static Vec::size_type GetArray( const std::string &ssKey,  _T ar[], Vec::size_type uSize)
	{
		DawgVar* pVar = GetVar(ssKey);
		if(pVar == NULL)
			return 0;
		return pVar->GetArray(ar, uSize);
	}
	template<class _T>
	static bool GetVector( const std::string &ssKey, std::vector<_T> &rVec)
	{
		DawgVar* pVar = GetVar(ssKey);
		if(pVar == NULL)
			return false;
		return pVar->GetVector(rVec);
	}

};

	//static Type GetType(double) { return tyNumber; }
	//static Type GetType(int)	{ return tyNumber; }
	//static Type GetType(bool)	{ return tyBool; }
	//static Type GetType(const Vec*) { return tyVector; }
	//static Type GetType(const std::string&) { return tyString; }
	//static Type GetType(const Tree*) {return tyTree; }

template <class Type> class SumValue
{
private:
	Type m_Sum;
public:
	SumValue () : m_Sum((Type)0) { }
	void operator ( ) ( const Type& elem ) {m_Sum += elem;}
    operator Type() const { return m_Sum; }
};

#endif
