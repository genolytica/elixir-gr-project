library(shiny)
library(shinyBS)
ui <- fluidPage(
    actionButton(
        inputId="showFastqc",
        label="Dataset FastQC",
        class="btn-sample-select"
    ),
    tags$head(tags$style(".modal-dialog{ width:1400px}")),
    bsModal(id="fqc", title="FastQC", trigger="click",
        includeHTML("www/GSE59017_multiqc_report.html")
    )  
)


server <- function(input, output, session) {
    observeEvent(input$showFastqc, {
        toggleModal(session, "fqc", toggle="toggle")
    })
}

shinyApp(ui, server)