process hifiasm {

    input:
    val (sample), val(ont), val(pb), val(hicr1), val(hicr2)

    output:
    val ("${sample}.h1"), path("*.fasta")

    script:

    def ontfiles = 
    def pbfiles = 

    """
    module load hifiasm
    


    """

}
