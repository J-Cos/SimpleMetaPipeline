    #if restart needed - specify where from
        RestartPipeline <- function(Pipeline, desiredlines) {
            lines <- scan(Pipeline, what=character(), sep="\n", skip=0, n=tail(desiredlines,1), quiet=TRUE)[desiredlines]
            tc <- textConnection(lines)
            source(tc)
            close(tc)
        }

        RestartPipeline(    Pipeline=file.path(path, "BioinformaticPipeline", "Pipeline", "CDLCLA_Pipeline.R"),
                            desiredlines=c(1:29, 123:180) )
