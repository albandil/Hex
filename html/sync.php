<?php

    // MIME type
    //   - According to http://vamdc.org/documents/standards/dataAccessProtocol/vamdctap.html
    //     the MIME type should be application/x-xsams+xml instead of simple text/xml.
    header('Content-Type: application/x-xsams+xml');

    // http://www.vamdc.eu/documents/standards/dataAccessProtocol/vamdctap.html
    //
    // Not to forget:
    // 
    //  1. Node must return valid XSAMS documents as defined by the latest VAMDC-XSAMS standard
    //     in any case when the response document is required.
    //  2. When database contains no data corresponding the query, node must respond with HTTP 204 Status Code,
    //     with empty response body, both for HEAD and GET requests to the /sync endpoint of VAMDC-TAP protocol.
    //  3. Node must respond to the vss2 query SELECT SPECIES with an XSAMS document containing only species
    //     information about all molecules, atoms and particles contained in node database, without states
    //     and processes data. Response should be given within a reasonable amount of time, no more than 30 seconds.
    //  4. Node must support gzip content encoding of transferred data, as defined in HTTP specification.
    //     This is employed to preserve the bandwidth and speed up the transfer of data, since XSAMS documents compress very well.
    //  5. Node must support HTTP HEAD requests to the TAP sync endpoint, giving sensible values in VAMDC-COUNT-* headers.
    //     Response to HEAD requests should be generated within a reasonable amount of time, no more than 30 seconds.
    //     Values may be inaccurate, but should give a view on how much data will be returned by the node for GET request.
    //  6. If node contains transitional data, it must support queries by RadTransWavelength, defining transition wavelength
    //     in vacuum in Angstroms. Use of this keyword is a common convention for clients querying the transition data.
    //  7. Node must provide sensible sample queries in Capabilities registration.

    // Time of the last modification of the database. Also good for caching.
    $time = time() - 60; // or filemtime($fn), etc
    header('Last-Modified: '.gmdate('D, d M Y H:i:s', $time).' GMT');

    // VAMDC headers informing on the content of the database.
    header('VAMDC-COUNT-ATOMS: 1');
    header('VAMDC-COUNT-MOLECULES: 0');
    header('VAMDC-COUNT-SPECIES: 2');
    header('VAMDC-COUNT-SOURCES: 1');
    header('VAMDC-COUNT-STATES: 100');
    header('VAMDC-COUNT-COLLISIONS: 1000');
    header('VAMDC-COUNT-RADIATIVE: 0');
    header('VAMDC-COUNT-NONRADIATIVE: 0');

    // Estimation of the size of the response document (MiB).
    header('VAMDC-APPROX-SIZE: 10');

    // Truncation of the response.
    // header('VAMDC-TRUNCATED: 2.9 %');

?>
<?xml version="1.0" encoding="UTF-8"?>
<XSAMSData xmlns="http://vamdc.org/xml/xsams/1.0"
           xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
           xsi:schemaLocation="http://vamdc.org/xml/xsams/1.0 http://vamdc.org/xml/xsams/1.0">
<?php
   
    error_reporting(-1);
    
    function hStateNL ($n, $l)
    {
        $symbols = array ( "s", "p", "d", "f", "g", "h", "i", "k", "l", "m", "o" );
        if ($l <= 10)
            return $n . $symbols[$l];
        else
            return $n . "(" . $l . ")";
    }
    
    function hStateNLM ($n, $l, $m)
    {
        return hStateNL($n,$l) . $m;
    }
    
    function IAEA_code ($ni, $li, $mi, $nf, $lf, $mf)
    {
        if ($ni == $nf and $li == $lf and $mi == $mf)
            return "EEL";
        if ($ni < $nf)
            return "EEX";
        if ($ni > $nf)
            return "EDX";
        return "EEX";
    }
    
    function collision_code ($ni, $li, $mi, $nf, $lf, $mf)
    {
        if ($ni == $nf and $li == $lf and $mi == $mf)
            return "elas";
        if ($ni < $nf)
            return "exci";
        if ($ni > $nf)
            return "deex";
        return "inel";
    }
    
    function expand_unqualified ($where_clause)
    {
        // (?<!\.) ... Negative lookbehind assertion : only match if '.' does not preceed the next pattern.
        $rex = "/(?<!\.)(StateEnergy|AtomStateTotalAngMom|AtomStateMagneticQuantumNumber|ParticleName|AtomSymbol)/i";
        
        while (preg_match($rex, $where_clause, $matches))
        {
            $reactant_part = preg_replace($rex, "reactant." . $matches[1], $where_clause);
            $product_part  = preg_replace($rex, "product." . $matches[1], $where_clause);
            $where_clause  = "((" . $reactant_part . ") AND (" . $product_part . "))";
        }
        return $where_clause;
    }
    
    function writeSpeciesSection ($db, $where_clause, $req_AtomStates, $req_Atoms, $req_Particles, $req_Species, $req_States)
    {
        $atom_string = Null;
        $particle_string = Null;
        
        if ($req_Species or $req_Atoms or $req_Particles or $req_AtomStates or $req_States)
        {
            if ($req_Species or $req_Atoms or $req_AtomStates or $req_States)
            {
                $statement1 = $db->prepare("SELECT DISTINCT reactant_AtomSymbol FROM vamdc_view " . $where_clause . ";");
                if (!$statement1)
                    return False;
                    
                $result1 = $statement1->execute();
                if (!$result1)
                    return False;
                
                $arr = $result1->fetchArray(SQLITE3_NUM);
                if ($arr)
                {
                    $atom_string =                "        <Atoms>\n";
                    $atom_string = $atom_string . "            <Atom>\n";
                    $atom_string = $atom_string . "                <ChemicalElement>\n";
                    $atom_string = $atom_string . "                    <NuclearCharge>1</NuclearCharge>\n";
                    $atom_string = $atom_string . "                    <ElementSymbol>H</ElementSymbol>\n";
                    $atom_string = $atom_string . "                </ChemicalElement>\n";
                    $atom_string = $atom_string . "                <Isotope>\n";
                    $atom_string = $atom_string . "                    <Ion speciesID = \"X.H\">\n";
                    $atom_string = $atom_string . "                        <IonCharge>0</IonCharge>\n";
                    
                    if ($req_AtomStates or $req_States)
                    {
                        $statement2 = $db->prepare("SELECT DISTINCT ni, li, mi FROM vamdc_view " . $where_clause . " UNION SELECT DISTINCT nf, lf, mf FROM vamdc_view " . $where_clause . ";");
                        if (!$statement2)
                            return False;
                        
                        $result2 = $statement2->execute();
                        if (!$result2)
                            return False;
                        
                        $arr = $result2->fetchArray(SQLITE3_NUM);
                        while ($arr)
                        {
                            $n = intval($arr[0]);
                            $l = intval($arr[1]);
                            $m = intval($arr[2]);
                            
                            $atom_string = $atom_string . "                        <AtomicState stateID = \"S." . hStateNLM($n,$l,$m) . "\">\n";
                            $atom_string = $atom_string . "                            <Description>" . hStateNL($n,$l) . " (m = " . $m . ")</Description>\n";
                            $atom_string = $atom_string . "                            <AtomicNumericalData>\n";
                            $atom_string = $atom_string . "                                <StateEnergy><Value units = \"Ry\">" . -1.0/($n*$n) . "</Value></StateEnergy>\n";
                            $atom_string = $atom_string . "                            </AtomicNumericalData>\n";
                            $atom_string = $atom_string . "                            <AtomicQuantumNumbers>\n";
                            $atom_string = $atom_string . "                                <TotalAngularMomentum>" . $l . "</TotalAngularMomentum>\n";
                            $atom_string = $atom_string . "                                <MagneticQuantumNumber>" . $m . "</MagneticQuantumNumber>\n";
                            $atom_string = $atom_string . "                            </AtomicQuantumNumbers>\n";
                            $atom_string = $atom_string . "                        </AtomicState>\n";
                            
                            $arr = $result2->fetchArray(SQLITE3_NUM);
                        }
                    }
                    
                    $atom_string = $atom_string . "                        <InChIKey>YZCKVEUIGOORGS-UHFFFAOYSA-N</InChIKey>\n";
                    $atom_string = $atom_string . "                    </Ion>\n";
                    $atom_string = $atom_string . "                </Isotope>\n";
                    $atom_string = $atom_string . "            </Atom>\n";
                    $atom_string = $atom_string . "        </Atoms>\n";
                }
            }
            
            if ($req_Species or $req_Particles)
            {
                $statement3 = $db->prepare("SELECT DISTINCT reactant_ParticleName FROM vamdc_view " . $where_clause . ";");
                if (!$statement3)
                    return False;
                
                $result3 = $statement3->execute();
                if (!$result3)
                    return False;
                
                $arr = $result3->fetchArray(SQLITE3_NUM);
                if ($arr)
                {
                    $particle_string =                    "        <Particles>\n";
                    $particle_string = $particle_string . "            <Particle speciesID = \"X.e\" name = \"electron\"/>\n";
                    $particle_string = $particle_string . "        </Particles>\n";
                }
            }
            
            if ($atom_string or $particle_string)
            {
                echo "    <Species>\n";
                echo $atom_string;
                echo $particle_string;
                echo "    </Species>\n";
                return True;
            }
        }
        
        return False;
    }
    
    function writeProcessesSection ($db, $where_clause, $req_Collisions, $req_Processes)
    {
        if ($req_Collisions or $req_Processes)
        {
            $statement = $db->prepare("SELECT DISTINCT ni, li, mi, nf, lf, mf FROM vamdc_view " . $where_clause . ";");
            if (!$statement)
                return False;
            
            $result = $statement->execute();
            if (!$result)
                return False;
            
            $substatement = $db->prepare("SELECT Ei, SUM(sigma), ni, li, mi, nf, lf, mf FROM ics WHERE ni = :ni AND li = :li AND mi = :mi AND nf = :nf AND lf = :lf AND mf = :mf GROUP BY ni,li,mi,nf,lf,mf,Ei;");
            if (!$substatement)
                return False;
            
            $substatement->bindParam(':ni', $ni, SQLITE3_INTEGER);
            $substatement->bindParam(':li', $li, SQLITE3_INTEGER);
            $substatement->bindParam(':mi', $mi, SQLITE3_INTEGER);
            $substatement->bindParam(':nf', $nf, SQLITE3_INTEGER);
            $substatement->bindParam(':lf', $lf, SQLITE3_INTEGER);
            $substatement->bindParam(':mf', $mf, SQLITE3_INTEGER);
            
            $arr = $result->fetchArray(SQLITE3_NUM);
            if ($arr)
            {
                $nodata = True;
                while ($arr)
                {
                    $ni = intval($arr[0]);
                    $li = intval($arr[1]);
                    $mi = intval($arr[2]);
                    $nf = intval($arr[3]);
                    $lf = intval($arr[4]);
                    $mf = intval($arr[5]);
                    
                    $substatement->reset();
                    $subresult = $substatement->execute();
                    if (!$subresult)
                        continue;
                    
                    $count = 0;
                    $subarr = $subresult->fetchArray(SQLITE3_NUM);
                    $energy_data = "";
                    $sigma_data = "";
                    while ($subarr)
                    {
                        $count++;
                        $energy_data = $energy_data . " " . $subarr[0];
                        $sigma_data = $sigma_data . " " . $subarr[1];
                        $subarr = $subresult->fetchArray(SQLITE3_NUM);
                    }
                    
                    if (!$count)
                        continue;
                    
                    if ($nodata)
                    {
                        echo "    <Processes>\n";
                        echo "        <Collisions>\n";
                        $nodata = False;
                    }
                    
                    echo "            <CollisionalTransition id = \"P." . hStateNLM($ni,$li,$mi) . "." . hStateNLM($nf,$lf,$mf) . "\">\n";
                    echo "                <ProcessClass>\n";
                    echo "                    <Code>" . collision_code($ni,$li,$mi,$nf,$lf,$mf) . "</Code>\n";
                    echo "                    <IAEACode>" . IAEA_code($ni,$li,$mi,$mf,$lf,$mf) . "</IAEACode>\n";
                    echo "                </ProcessClass>\n";
                    echo "                <Reactant>\n";
                    echo "                    <SpeciesRef>X.H</SpeciesRef>\n";
                    echo "                    <StateRef>S." . hStateNLM($ni,$li,$mi) . "</StateRef>\n";
                    echo "                </Reactant>\n";
                    echo "                <Reactant>\n";
                    echo "                    <SpeciesRef>X.e</SpeciesRef>\n";
                    echo "                </Reactant>\n";
                    
                    if ($nf == 0)
                    {
                        echo "                <Product>\n";
                        echo "                    <SpeciesRef>X.p</SpeciesRef>\n";
                        echo "                </Product>\n";
                        echo "                <Product>\n";
                        echo "                    <SpeciesRef>X.e</SpeciesRef>\n";
                        echo "                </Product>\n";
                        
                        $threshold = 1.0 / (ni * ni);
                    }
                    else
                    {
                        echo "                <Product>\n";
                        echo "                    <SpeciesRef>X.H</SpeciesRef>\n";
                        echo "                    <StateRef>S." . hStateNLM($nf,$lf,$mf) . "</StateRef>\n";
                        echo "                </Product>\n";
                        echo "                <Product>\n";
                        echo "                    <SpeciesRef>X.e</SpeciesRef>\n";
                        echo "                </Product>\n";
                        
                        if ($ni < $nf)
                            $threshold = 0.0;
                        else
                            $threshold = 1.0 / ($ni * $ni) - 1.0 / ($nf * $nf);
                    }
                    
                    echo "                <Threshold><Value units = \"Ry\">" . $threshold . "</Value></Threshold>\n";
                    echo "                <DataSets>\n";
                    echo "                    <DataSet dataDescription = \"crossSection\">\n";
                    echo "                        <TabulatedData>\n";
                    echo "                            <Description>Scattering cross section for reaction e⁻ + H(" . $ni . "," . $li . "," . $mi . ") -> e⁻ + H(" . $nf . "," . $lf . "," . $mf . ")</Description>\n";
                    echo "                            <X units = \"Ry\"><DataList count = \"$count\">" . $energy_data . "</DataList></X>\n";
                    echo "                            <Y units = \"au\"><DataList>" . $sigma_data . "</DataList></Y>\n";
                    echo "                            <ReferenceFrame>TargetFrame</ReferenceFrame>\n";
                    echo "                        </TabulatedData>\n";
                    echo "                    </DataSet>\n";
                    echo "                </DataSets>\n";
                    echo "            </CollisionalTransition>\n";
                    
                    $arr = $result->fetchArray(SQLITE3_NUM);
                }
                
                if (!$nodata)
                {
                    echo "        </Collisions>\n";
                    echo "    </Processes>\n";
                }
            }
        }
        
        return False;
    }
    
    function writeSourcesSection ($db, $where_clause, $req_Sources)
    {
//         if ($req_Sources)
        {
            $statement = $db->prepare("SELECT DISTINCT ni,li,mi,nf,lf,mf FROM vamdc_view " . $where_clause . ";");
            if (!$statement)
                return False;
            
            $result = $statement->execute();
            if (!$result)
                return False;
            
//             if ($result->fetchArray(SQLITE3_NUM))
            {
                echo "    <Sources>\n";
                echo "        <Source sourceID=\"B.HexDataI\">\n";
                echo "            <Category>journal</Category>\n";
                echo "            <SourceName>Comput. Phys. Commun. 213</SourceName>\n";
                echo "            <Year>2017</Year>\n";
                echo "            <Authors>\n";
                echo "                <Author><Name>Benda J.</Name></Author>\n";
                echo "                <Author><Name>Houfek K.</Name></Author>\n";
                echo "            </Authors>\n";
                echo "        </Source>\n";
                echo "    </Sources>\n";
                
                return True;
            }
        }
        
        return False;
    }
    
    function writeMethodsSection ($db, $where_clause, $req_Methods)
    {
//         if ($req_Methods)
        {
            $statement = $db->prepare("SELECT DISTINCT ni,li,mi,nf,lf,mf FROM vamdc_view " . $where_clause . ";");
            if (!$statement)
                return False;
            
            $result = $statement->execute();
            if (!$result)
                return False;
            
//             if ($result->fetchArray(SQLITE3_NUM))
            {
                echo "    <Methods>\n";
                echo "        <Method methodID=\"M.calc\">\n";
                echo "            <Category>theory</Category>\n";
                echo "            <Description>Solution of Schroedinger equation in B-spline basis with exterior complex scaling boundary condition.</Description>\n";
                echo "        </Method>\n";
                echo "    </Methods>\n";
                
                return True;
            }
        }
        
        return False;
    }

// -------------------------------------------------------------------------------------
// Analyze GET input
//

    // requested Requestables
    $req_AtomStates = False;
    $req_Atoms = False;
    $req_Collisions = False;
    $req_Functions = False;                 // ignored
    $req_Methods = False;
    $req_MoleculeBasisStates = False;       // ignored
    $req_MoleculeQuantumNumbers = False;    // ignored
    $req_MoleculeStates = False;            // ignored
    $req_Molecules = False;                 // ignored
    $req_NonRadiativeTransitions = False;   // ignored
    $req_Particles = False;
    $req_Processes = False;
    $req_RadiativeCrossSections = False;    // ignored
    $req_RadiativeTransitions = False;      // ignored
    $req_Solids = False;                    // ignored
    $req_Sources = False;
    $req_Species = False;
    $req_States = False;
    
    echo "    <!-- REQUEST: " . $_GET["REQUEST"] . " -->\n";
    echo "    <!-- LANG: " . $_GET["LANG"] . " -->\n";
    echo "    <!-- FORMAT: " . $_GET["FORMAT"] . " -->\n";
    echo "    <!-- QUERY: " . $_GET["QUERY"] . " -->\n";
    
    // Allowed REQUESTs are : doQuery.
    if ($_GET["REQUEST"] != "doQuery")
    {
        http_response_code(400); // Bad request with malformed query string or missing restrictable.
        exit;
    }

    // Allowed LANGs are : VSS2.
    if (strtoupper($_GET["LANG"]) != "VSS2")
    {
        http_response_code(400); // Bad request with malformed query string or missing restrictable.
        exit;
    }

    // Allowed FORMATs are : XSAMS.
    if (strtoupper($_GET["FORMAT"]) != "XSAMS")
    {
        http_response_code(400); // Bad request with malformed query string or missing restrictable.
        exit;
    }
    
    // All queries are allowed.

// -------------------------------------------------------------------------------------
// Parse the query
//

    // convert query to lower case
    $query = strtolower($_GET["QUERY"]);
    
    // get first token of the query : must be "select"
    $first_token = strtok($query, " ");
    if ($first_token != "select")
    {
        http_response_code(400); // Bad request with malformed query string or missing restrictable.
        exit;
    }
    
    // read further tokens, until "where" is found
    do
    {
        $next_token = strtok(" ,");
        switch ($next_token)
        {
            case "atomstates":              $req_AtomStates = 1;                break;
            case "atoms":                   $req_Atoms = 1;                     break;
            case "collisions":              $req_Collisions = 1;                break;
            case "functions":               $req_Functions = 1;                 break;
            case "methods":                 $req_Methods = 1;                   break;
            case "moleculebasisstates":     $req_MoleculeBasisStates = 1;       break;
            case "moleculequantumnumbers":  $req_MoleculeQuantumNumbers = 1;    break;
            case "moleculestates":          $req_MoleculeStates = 1;            break;
            case "molecules":               $req_Molecules = 1;                 break;
            case "nonradiativetransitions": $req_NonRadiativeTransitions = 1;   break;
            case "particles":               $req_Particles = 1;                 break;
            case "processes":               $req_Processes = 1;                 break;
            case "radiativecrosssections":  $req_RadiativeCrossSections = 1;    break;
            case "radiativetransitions":    $req_RadiativeTransitions = 1;      break;
            case "solids":                  $req_Solids = 1;                    break;
            case "sources":                 $req_Sources = 1;                   break;
            case "species":                 $req_Species = 1;                   break;
            case "states":                  $req_States = 1;                    break;
            case "*":
            case "all":
                $req_AtomStates = 1;
                $req_Atoms = 1;
                $req_Collisions = 1;
                $req_Functions = 1;
                $req_Methods = 1;
                $req_MoleculeBasisStates = 1;
                $req_MoleculeQuantumNumbers = 1;
                $req_MoleculeStates = 1;
                $req_Molecules = 1;
                $req_NonRadiativeTransitions = 1;
                $req_Particles = 1;
                $req_Processes = 1;
                $req_RadiativeCrossSections = 1;
                $req_RadiativeTransitions = 1;
                $req_Solids = 1;
                $req_Sources = 1;
                $req_Species = 1;
                $req_States = 1;
                break;
            case "where":
                break;
            default:
                if ($next_token != "")
                {
                    http_response_code(400); // Bad request with malformed query string or missing restrictable.
                    exit;
                }
        }
    }
    while ($next_token != null and $next_token != "where");
    
    // Find position of "where" in the query.
    $where_pos = strpos($query, "where");
    
    // If found, process the where clause.
    $where_clause = "";
    if ($where_pos)
    {
        // Get the whole "where" clause.
        $where_clause = substr($query, $where_pos + 5);
        
        // Expand unqualified Restrictables to both sides ("reactant_" & "product_")
        $where_clause = expand_unqualified($where_clause);
        
        // Translate the Restrictables:
        //      'target\.'              -> 'reactant_'
        $where_clause = preg_replace("/target\\./", "reactant_", $where_clause);
        //      'collider\.'            -> 'reactant_'
        $where_clause = preg_replace("/collider\\./", "reactant_", $where_clause);
        //      'reactant[0-9a-z]+\.'   -> 'reactant_'
        $where_clause = preg_replace("/reactant[0-9a-z]*\\./", "reactant_", $where_clause);
        //      'product[0-9a-z]+\.'    -> 'product_'
        $where_clause = preg_replace("/product[0-9a-z]*\\./", "product_", $where_clause);
        //      'lower\.'               -> 'lower_'
        $where_clause = preg_replace("/lower\\./", "lower_", $where_clause);
        //      'upper\.'               -> 'upper_'
        $where_clause = preg_replace("/upper\\./", "upper_", $where_clause);
        
        // Add 'where' keyword.
        $where_clause = "where " . $where_clause;
    }
    
    echo "    <!-- Adapted WHERE clause: " . $where_clause . " -->\n";
    
// -------------------------------------------------------------------------------------
// Process the query
//

    // Open the Hex database.
    include "paths.inc";
    $db = new SQLite3 ($hexdbdat, SQLITE3_OPEN_READONLY);
    
    // Write all sections of the XSAMS file.
    $valid_species_data = writeSpeciesSection($db, $where_clause, $req_AtomStates, $req_Atoms, $req_Particles, $req_Species, $req_States);
    $valid_process_data = writeProcessesSection($db, $where_clause, $req_Collisions, $req_Processes);
    $valid_sources_data = writeSourcesSection($db, $where_clause, $req_Sources);
    $valid_methods_data = writeMethodsSection($db, $where_clause, $req_Methods);
    
    // Close the database.
    $db->close();
    
    // Valid data case:
    if ($valid_species_data or $valid_process_data or $valid_sources_data or $valid_methods_data)
        http_response_code(200);
    
    // Request processed, but no matching data found.
    else
        http_response_code(204);

// Bad request with malformed query string or missing restrictable.
// http_response_code(400);

// Internal crash.
// http_response_code(500);

?>
</XSAMSData>
