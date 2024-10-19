pcf_function <- function(hf, i, j){
    #
    # hf is a hyperframe with a column called 'coord' containing the ppp
    # i mark 1
    # j mark 2
    #
    
    # Create a safe function in case there are errors in the hyperframe
    safe_envelope = possibly(.f = envelope, otherwise=NULL)
    
    # Apply inhomogeneous PCF with Markov Chain envelope
    pcf <- map(hf$coord, safe_envelope, pcfcross.inhom, i=i, j=j, divisor = "d", correction = "isotropic", r=seq(0,30,0.2))
    
    # Creates a tibble with the data from PCF
    pcf_tibble <- map(pcf, tibble) %>% 
        enframe(name = "rep") %>%
        unnest(cols = "value")
    
    return(pcf_tibble)
}