import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const workbookPath = process.argv[2];
const sheetName = process.argv[3];
const range = process.argv[4];

if (!workbookPath) {
  throw new Error("usage: inspect_supplementary.mjs <workbook.xlsx> [sheet] [range]");
}

const input = await FileBlob.load(workbookPath);
const workbook = await SpreadsheetFile.importXlsx(input);

if (!sheetName) {
  const result = await workbook.inspect({
    kind: "workbook,sheet,table",
    include: "id,name",
    maxChars: 12000,
    tableMaxRows: 8,
    tableMaxCols: 12,
    tableMaxCellChars: 120,
  });
  process.stdout.write(`${result.ndjson}\n`);
} else {
  const result = await workbook.inspect({
    kind: "region",
    sheetId: sheetName,
    range: range ?? "A1:Z80",
    include: "values,formulas",
    maxChars: 30000,
    tableMaxRows: 100,
    tableMaxCols: 30,
    tableMaxCellChars: 200,
  });
  process.stdout.write(`${result.ndjson}\n`);
}
