initial <- list()

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "Gm37180", Search = "", SearchColumns = c("",
"", "", "", "", "", ""), HiddenColumns = character(0), VersionInfo = list(
    iSEE = structure(list(c(2L, 12L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(RowDataTable = 1L), PanelHeight = 600L,
    PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data",
    XAxisColumnData = "Condition", XAxisFeatureName = "4933401J01Rik",
    XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
    YAxisFeatureName = "Gm37180", YAxisFeatureSource = "RowDataTable1",
    YAxisFeatureDynamicSource = TRUE, FacetRowByColData = "Sex",
    FacetColumnByColData = "Sex", ColorByColumnData = "SAMPLE_ID",
    ColorByFeatureNameAssay = "counts", ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "Condition", SizeByColumnData = "ERCCsumLogErr",
    TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None",
    ColorBy = "Column data", ColorByDefaultColor = "#000000",
    ColorByFeatureName = "4933401J01Rik", ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "Sample.10_1",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "Column data", SizeBy = "None", SelectionAlpha = 0.1,
    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
    VisualChoices = c("Color", "Size", "Text"), ContourAdd = FALSE,
    ContourColor = "#0000FF", PointSize = 5, PointAlpha = 1,
    Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
    CustomLabelsText = "Sample.10_1", FontSize = 1, LegendPointSize = 1,
    LegendPosition = "Bottom", HoverInfo = TRUE, LabelCenters = FALSE,
    LabelCentersBy = "Sex", LabelCentersColor = "#000000", VersionInfo = list(
        iSEE = structure(list(c(2L, 12L, 0L)), class = c("package_version",
        "numeric_version"))), PanelId = c(FeatureAssayPlot = 1L),
    PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list())
