-- .headers on
.mode tabs

-- Get all available initial states.
CREATE TEMP TABLE States AS 
	SELECT DISTINCT
		ni, li, mi
	FROM
		tmat
	ORDER BY
		ni, li, mi;

-- For every initial state and every incident energy, get maximal partial wave.
CREATE TEMP TABLE MaxLs AS
	SELECT
		S.ni AS ni, S.li AS li, S.mi AS mi, T.Ei AS Ei, MAX(T.L) AS MaxL
	FROM
		States AS S, tmat AS T 
	WHERE 
		S.ni = T.ni AND S.li = T.li AND S.mi = T.mi
	GROUP BY
		S.ni, S.li, S.mi, T.Ei;

-- Collapse energy intervals of constant MaxL
CREATE TEMP TABLE StackedLs AS
	SELECT 
		LThis.ni, LThis.li, LThis.mi, 
		LThis.Ei AS ThisE, LThis.MaxL AS ThisL
	FROM 
		MaxLs AS LThis 
		LEFT OUTER JOIN MaxLs AS LPrev
		ON LThis.rowid = LPrev.rowid + 1 
			AND LThis.ni = LPrev.ni 
			AND LThis.li = LPrev.li 
			AND LThis.mi = LPrev.mi
		LEFT OUTER JOIN MaxLs AS LNext
		ON LThis.rowid = LNext.rowid - 1 
			AND LThis.ni = LNext.ni 
			AND LThis.li = LNext.li 
			AND LThis.mi = LNext.mi
	WHERE 
		LNext.MaxL IS NULL OR LPrev.MaxL IS NULL OR LPrev.MaxL <> LNext.MaxL
	ORDER BY 
		LThis.ni, LThis.li, LThis.mi, ThisE;

-- Print out the first and lase energies for every constant interval.
SELECT 
	T1.ni, T1.li, T1.mi, T1.ThisE AS MinE, (CASE WHEN T2.ThisE IS NULL THEN T1.ThisE ELSE T2.ThisE END) AS MaxE, T1.ThisL AS MaxL
FROM 
	StackedLs AS T1
	LEFT OUTER JOIN StackedLs AS T2
	ON T1.rowid = T2.rowid - 1
		AND T1.ni = T2.ni
		AND T1.li = T2.li
		AND T1.mi = T2.mi
	WHERE (T1.ThisL = T2.ThisL OR T2.ThisL IS NULL) AND MinE != MaxE;
