library(shiny)
#library(shinyBS)

ui <- fluidPage(
	tags$head(tags$style(".modal-dialog{ width:1400px}")),
	navbarPage(
        id="test",
        title="MultiQC test",
        tabPanel("Button",
            actionButton(
				inputId="showFastqc",
				label="Dataset FastQC",
				class="btn-sample-select"
			)
        ),
        tabPanel("Static"#,
            #includeHTML("www/GSE59017_multiqc_report.html")
        )
    )
)

server <- function(input, output, session) {
	observeEvent(input$showFastqc, {
		showModal(modalDialog(
			title="Test",
			#includeHTML("www/GSE59017_multiqc_report.html")
			tags$iframe(
				src="GSE59017_multiqc_report.html",
				width="100%",height=640,
				#sandbox=paste("allow-forms allow-popups allow-pointer-lock",
				#	"allow-same-origin allow-scripts allow-top-navigation"),
				frameborder=0
			),
			size="l"
		))
	})
}

shinyApp(ui, server)
