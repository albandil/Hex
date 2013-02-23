-- .headers on
.mode tabs

-- CREATE TABLE tmat (
-- 	ni INTEGER,
-- 	li INTEGER,
-- 	mi INTEGER,
-- 	Ei DOUBLE PRECISION,
-- 	L INTEGER
-- );
-- 
-- INSERT INTO tmat VALUES (1, 0, 0, 0.1, 0);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.1, 1);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.1, 2);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.1, 3);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.2, 0);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.2, 1);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.2, 2);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.2, 3);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.3, 0);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.3, 1);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.3, 2);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.3, 3);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.4, 0);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.4, 1);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.4, 2);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.5, 0);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.5, 1);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.5, 2);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.6, 0);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.6, 1);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.6, 2);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.6, 3);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.7, 0);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.7, 1);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.7, 2);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.7, 3);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.8, 0);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.8, 1);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.8, 2);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.8, 3);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.9, 0);
-- INSERT INTO tmat VALUES (1, 0, 0, 0.9, 1);
-- INSERT INTO tmat VALUES (2, 0, 0, 0.5, 0);
-- INSERT INTO tmat VALUES (2, 0, 0, 0.5, 1);

CREATE TEMP TABLE MaxLs AS
	SELECT 
		ni, li, mi, Ei, MAX(L) AS MaxL
	FROM
		tmat
	GROUP BY
		ni,li,mi,Ei
	ORDER BY
		Ei ASC;

-- SELECT * FROM MaxLs;

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

-- SELECT * FROM StackedLs;

CREATE TEMP TABLE SStackedLs AS
	SELECT 
		T1.ni, T1.li, T1.mi, T1.ThisE AS MinE, (CASE WHEN T2.ThisE IS NULL THEN T1.ThisE ELSE T2.ThisE END) AS MaxE, T1.ThisL AS MaxL
	FROM 
		StackedLs AS T1
		LEFT OUTER JOIN StackedLs AS T2
		ON T1.rowid = T2.rowid - 1
			AND T1.ni = T2.ni
			AND T1.li = T2.li
			AND T1.mi = T2.mi
		WHERE T1.ThisL = T2.ThisL OR T2.ThisL IS NULL
;

SELECT * FROM SStackedLs;
