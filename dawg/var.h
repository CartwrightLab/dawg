#ifndef DAWG_VAR_H
#define DAWG_VAR_H

#include "dawg.h"
#include "tree.h"

class MapSsToVar;

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
	explicit DawgVar(NewickNode* p) : m_tyType(tyTree), m_ptrData(p) { }
	NewickNode*	GetTree() const { return m_ptrData; }
	bool		Get(NewickNode *&rTree) const;
	void		Set(NewickNode *pTree);

protected:
	Type m_tyType;
	union
	{
		double m_dData;
		bool   m_bData;
		const Vec*   m_pvData;
		const std::string* m_pssData;
		NewickNode*  m_ptrData;
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

class MapSsToVar : public std::map<std::string, DawgVar*>
{
public:	
	virtual ~MapSsToVar()
	{
		for(iterator it = begin(); it != end(); ++it)
			if(it->second)
				delete it->second;
		clear();
	}
};

#endif //DAWG_VAR_H