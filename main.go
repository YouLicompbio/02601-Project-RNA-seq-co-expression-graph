package main

import (
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
	"flag"
)
/*example input:
const (
	// input1: GCT data path
	gctDataFile = "gene_reads_v10_thyroid.gct.gz"

	// input2: GTF annotation path (length of genes, for TPM)
	gtfAnnotationFile = "gencode.v44.annotation.gtf.gz"

	// output: matrix cleaned
	outputMatrixFile = "clean_thyroid_matrix.csv"

	// Filtering parameter
	// We delete a gene if it's expressed in 90% samples (based on log2(TPM+1) < 1)
	lowExpressionThreshold = 0.9

	// delete genes of low variance percentile
	lowVariancePercentile = 0.25

	//soft threshold: beta
	softPowerBeta = 6.0
)
*/

type Config struct {
	GCTPath       string
	GTFPath       string
	OutputPrefix  string
	LowExpThresh  float64
	LowVarPercent float64
	SoftPowerBeta float64
}

func init() {
	// 定义命令行参数，RShiny 将通过这些参数控制 Go 程序
	flag.StringVar(&config.GCTPath, "gct", "gene_reads.gct.gz", "Path to the input GCT file (gene counts)")
	flag.StringVar(&config.GTFPath, "gtf", "annotation.gtf.gz", "Path to the GTF annotation file")
	flag.StringVar(&config.OutputPrefix, "out", "output", "Prefix for the output files (e.g., 'results/run1')")
	
	flag.Float64Var(&config.LowExpThresh, "min_exp", 0.9, "Threshold for low expression filtering (log2CPM)")
	flag.Float64Var(&config.LowVarPercent, "min_var", 0.25, "Percentile for low variance filtering (0-1)")
	flag.Float64Var(&config.SoftPowerBeta, "beta", 6.0, "Soft thresholding power (beta) for network construction")
	
	// 自定义 Usage 信息，提升文档质量
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage of %s:\n", os.Args[0])
		fmt.Fprintf(os.Stderr, "This program constructs a gene co-expression network from RNA-seq data.\n")
		fmt.Fprintf(os.Stderr, "\nParameters:\n")
		flag.PrintDefaults()
	}
}

var config Config


func main() {
	flag.Parse() // 解析参数

	// ---------------------------------------------------------
	// PHASE 1: Preprocessing (Gene x Sample)
	// ---------------------------------------------------------
	log.Println("Phase 1: Preprocessing the data...")
	
	log.Println("Parsing GTF annotation...")
	geneLengthsKB, err := parseGTFToLengths(config.GTFPath)
	if err != nil {
		log.Fatalf("Failed to parse GTF: %v", err)
	}

	log.Println("Preprocessing GCT raw counts...")
	finalMatrix, finalGeneList, finalSampleList, err := processGCTFile(
		config.GCTPath,
		geneLengthsKB,
		config.LowExpThresh,
		config.LowVarPercent,
	)
	if err != nil {
		log.Fatalf("Failed to process GCT: %v", err)
	}

	// [Phase 1Output] 这里还是用 writeOutputCSV，因为它是 Gene x Sample
	cleanMatrixFile := config.OutputPrefix + "_clean_matrix.csv"
	log.Println("Saving clean matrix to:", cleanMatrixFile)
	err = writeOutputCSV(cleanMatrixFile, finalMatrix, finalGeneList, finalSampleList)
	if err != nil {
		log.Fatalf("Failed writing clean matrix: %v", err)
	}

	// ---------------------------------------------------------
	// PHASE 2: Correlation (Gene x Gene)
	// ---------------------------------------------------------
	log.Println("Phase 2: Calculating Correlation Matrix...")
	correlationMatrix, err := RunPhase2(finalMatrix, finalGeneList)
	if err != nil {
		log.Fatalf("Failed Phase 2: %v", err)
	}

	// [使用你的函数] 保存 Correlation
	corrFile := config.OutputPrefix + "_correlation.csv"
	log.Println("Saving Correlation Matrix...")
	err = writeCorrelationMatrix(corrFile, correlationMatrix, finalGeneList)
	if err != nil {
		log.Printf("warning: failed to save correlation: %v", err)
	}

	// ---------------------------------------------------------
	// PHASE 3: Adjacency (Gene x Gene)
	// ---------------------------------------------------------
	log.Printf("Phase 3: Calculating Adjacency (Beta=%.1f)...", config.SoftPowerBeta)
	adjacencyMatrix := CalculateAdjacencyMatrix(correlationMatrix, config.SoftPowerBeta)
	
	// [使用你的函数] 保存 Adjacency
	adjFile := config.OutputPrefix + "_adjacency.csv"
	log.Println("Saving Adjacency Matrix...")
	err = writeCorrelationMatrix(adjFile, adjacencyMatrix, finalGeneList)
	if err != nil {
		log.Printf("warning: failed to save adjacency: %v", err)
	}

	// ---------------------------------------------------------
	// PHASE 4: TOM (Gene x Gene)
	// ---------------------------------------------------------
	log.Println("Phase 4: Calculating TOM...")
	tomMatrix := CalculateTOM(adjacencyMatrix)
	
	// [使用你的函数] 保存 TOM
	tomFile := config.OutputPrefix + "_tom_matrix.csv"
	log.Println("Saving TOM Matrix...")
	err = writeCorrelationMatrix(tomFile, tomMatrix, finalGeneList)
	if err != nil {
		log.Printf("warning: failed to save TOM: %v", err)
	}

	// ---------------------------------------------------------
	// PHASE 5: Dissimilarity (Gene x Gene) -> For RShiny
	// ---------------------------------------------------------
	log.Println("Phase 5: Calculating Dissimilarity (1-TOM)...")
	distMatrix := CalculateDissimilarity(tomMatrix)
	
	// [使用你的函数] 保存 Dissimilarity
	finalFile := config.OutputPrefix + "_dissimilarity.csv"
	log.Println("Saving Dissimilarity Matrix for clustering...")
	err = writeCorrelationMatrix(finalFile, distMatrix, finalGeneList)
	if err != nil {
		log.Fatalf("Failed to save final dissimilarity matrix: %v", err)
	}

	log.Println("DONE! Pipeline finished.")
}


// writeOutputCSV takes in filepath, processed matrix, genelist, samplelist input
// It saves the matrix in local environment.
func writeOutputCSV(filePath string, matrix [][]float64, geneList []string, sampleList []string) error {
	file, err := os.Create(filePath)
	if err != nil {
		return err
	}
	defer file.Close()

	// Write in sample id
	// add a gene_id column
	header := append([]string{"gene_id"}, sampleList...)
	_, err = fmt.Fprintln(file, strings.Join(header, ","))
	if err != nil {
		return err
	}

	// write in the gene data
	for i, geneID := range geneList {
		row := matrix[i]
		rowStr := make([]string, len(row)+1)
		rowStr[0] = geneID
		for j, val := range row {
			rowStr[j+1] = strconv.FormatFloat(val, 'f', 6, 64) 
		}
		_, err = fmt.Fprintln(file, strings.Join(rowStr, ","))
		if err != nil {
			
			log.Printf("Failed to write %s, %v", geneID, err)
		}
	}
	return nil
}


