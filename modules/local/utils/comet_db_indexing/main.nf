process COMET_DB_INDEXING {
    tag "$database.baseName"
    label 'process_medium'
    label 'openms'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/bigbio/openms-tools-thirdparty-sif:2025.04.14' :
        'ghcr.io/bigbio/openms-tools-thirdparty:2025.04.14' }"

    input:
    tuple val(meta), path(database)

    output:
    path("${database}.idx"), emit: idx_file
    path "versions.yml",     emit: versions
    path "*.log",            emit: log

    script:
    def args = task.ext.args ?: ''

    // Handle fragment mass tolerance unit conversion (same logic as COMET module)
    if (meta.fragmentmasstoleranceunit == "ppm") {
        // Note: This uses an arbitrary rule to decide if it was hi-res or low-res
        // and uses Comet's defaults for bin size, in case unsupported unit "ppm" was given.
        if (meta.fragmentmasstolerance.toDouble() < 50) {
            bin_tol = 0.015
            bin_offset = 0.0
            inst = params.instrument ?: "high_res"
        } else {
            bin_tol = 0.50025
            bin_offset = 0.4
            inst = params.instrument ?: "low_res"
        }
        log.warn "The chosen search engine Comet does not support ppm fragment tolerances. We guessed a " + inst +
            " instrument and set the fragment_bin_tolerance to " + bin_tol
    } else {
        bin_tol = meta.fragmentmasstolerance.toDouble()
        bin_offset = bin_tol <= 0.05 ? 0.0 : 0.4
        if (!params.instrument)
        {
            inst = bin_tol <= 0.05 ? "high_res" : "low_res"
        } else {
            inst = params.instrument
        }
    }

    // Handle isotope error range
    def isoSlashComet = "0/1"
    if (params.isotope_error_range) {
        def isoRangeComet = params.isotope_error_range.split(",")
        isoSlashComet = ""
        for (c in isoRangeComet[0].toInteger()..isoRangeComet[1].toInteger()-1) {
            isoSlashComet += c + "/"
        }
        isoSlashComet += isoRangeComet[1]
    }

    // Handle enzyme compatibility with MSGF (same logic as COMET module)
    enzyme = meta.enzyme
    if (params.search_engines.contains("msgf")){
        if (meta.enzyme == "Trypsin") enzyme = "Trypsin/P"
        else if (meta.enzyme == "Arg-C") enzyme = "Arg-C/P"
        else if (meta.enzyme == "Asp-N") enzyme = "Arg-N/B"
        else if (meta.enzyme == "Chymotrypsin") enzyme = "Chymotrypsin/P"
        else if (meta.enzyme == "Lys-C") enzyme = "Lys-C/P"
    }

    // Handle enzyme termini
    num_enzyme_termini = ""
    if (meta.enzyme == "unspecific cleavage")
    {
        num_enzyme_termini = "none"
    }
    else if (params.num_enzyme_termini == "fully")
    {
        num_enzyme_termini = "full"
    }

    il_equiv = params.IL_equivalent ? "-PeptideIndexing:IL_equivalent" : ""
    met_excision = params.met_excision ? "-clip_nterm_methionine true" : ""

    """
    CometAdapter \\
        -database "${database}" \\
        -fragindex \\
        -instrument ${inst} \\
        -missed_cleavages $params.allowed_missed_cleavages \\
        -min_peptide_length $params.min_peptide_length \\
        -max_peptide_length $params.max_peptide_length \\
        -num_hits $params.num_hits \\
        -num_enzyme_termini $params.num_enzyme_termini \\
        -enzyme "${enzyme}" \\
        -isotope_error ${isoSlashComet} \\
        -precursor_charge $params.min_precursor_charge:$params.max_precursor_charge \\
        -fixed_modifications ${meta.fixedmodifications.tokenize(',').collect { "'$it'" }.join(" ") } \\
        -variable_modifications ${meta.variablemodifications.tokenize(',').collect { "'$it'" }.join(" ") } \\
        -max_variable_mods_in_peptide $params.max_mods \\
        -precursor_mass_tolerance $meta.precursormasstolerance \\
        -precursor_error_units $meta.precursormasstoleranceunit \\
        -fragment_mass_tolerance ${bin_tol} \\
        -fragment_bin_offset ${bin_offset} \\
        -minimum_peaks $params.min_peaks \\
        ${met_excision} \\
        ${il_equiv} \\
        -debug $params.db_debug \\
        -force \\
        $args \\
        2>&1 | tee ${database.baseName}_comet_idx.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CometAdapter: \$(CometAdapter 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
        Comet: \$(comet 2>&1 | grep -E "Comet version.*" | sed 's/ Comet version //g' | sed 's/"//g')
    END_VERSIONS
    """
}
