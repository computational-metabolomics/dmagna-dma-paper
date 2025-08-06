
x <- read.csv('output_original/daphnia_annotation_summary.csv')
y <- read.csv('output/daphnia_annotation_summary.csv')
compare(x, y)

x <- read.csv('output_original/daphnia_annotations_all_class_count.csv')
y <- read.csv('output/daphnia_annotations_all_class_count.csv')
compare(x, y)

x <- read.csv('output_original/daphnia_annotations_all_subclass_count.csv')
y <- read.csv('output/daphnia_annotations_all_subclass_count.csv')
compare(x, y)

x <- read.csv('output_original/daphnia_annotations_all_superclass_count.csv')
y <- read.csv('output/daphnia_annotations_all_superclass_count.csv')
compare(x, y)
