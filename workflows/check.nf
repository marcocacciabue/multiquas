/* 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Check Files and version 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow CHECK {

	// Require version

		if (!nextflow.version.matches('>=25.0.0')) {
			error ("ERROR: This workflow requires Nextflow version '25.0.0' or later. You are running '${nextflow.version}'.")
		}


	// Check FASTA inputs

		def fasta_extensions = ['fa', 'fna', 'fasta']
		def query = file("${params.input_ref}")
  
		if (query.exists()) {
			def file_extension = query.extension
				if (!fasta_extensions.contains(file_extension)) {
					error ("ERROR: File '${query}' does not have a '.fa', '.fna', or '.fasta' extension, please provide a FASTA file.")
				}
		} else {
			error ("ERROR: Reference file not found at '${query}'. Please check you have indicated the filepath correctly.")
		}

}
