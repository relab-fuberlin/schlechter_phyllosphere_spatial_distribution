
##AREA UNDER THE CURVE

library(DescTools)

# Load data
fractions <- read.csv(here('results', 'stat_K_fractions.csv'), header = TRUE, row.names = "X") %>% 
    tibble %>% 
    separate(col='syncom', into=c('synID', 'synC'), sep = '\\.', remove = FALSE)
fractions %>% str
fractions %>% head

auc_fractions <- fractions %>% 
    group_by(syncom, synID, dpi, type, strain) %>% 
    #group_by(synID, dpi, type, strain) %>% 
    summarise(auc = AUC(r, fraction, absolutearea = TRUE),
              .groups = 'drop')

auc_C <- auc_fractions %>% 
    filter(synID == 'C') %>% 
    #select(-synID) %>% 
    select(-syncom:-synID)

auc_fold_change <- auc_fractions %>% 
    filter(synID != 'C') %>% 
    inner_join(., auc_C, by = c('dpi', 'type', 'strain'), suffix = c('.inter', '.intra')) %>% 
    mutate(fold = auc.inter/auc.intra,
           log2FC = log2(fold))

auc_fold_change %>% 
    ggplot(aes(syncom, strain))+
    facet_grid(dpi~type)+
    geom_tile(colour='black', fill='white')+
    geom_point(aes(fill=log2FC), shape = 21, size = 12)+
    scale_fill_gradientn(colours = wes_palette("Zissou1"))+
    coord_fixed()

auc_fold_change %>% 
    ggplot(aes(dpi, log2FC))+
    facet_grid(type~strain)+
    geom_point()

auc_fractions %>% 
    ggplot(aes(auc))+
    geom_histogram()

lmauc = lm(auc ~ synID + dpi + type + strain, data = auc_fractions)
shapiro.test(rstandard(lmauc))
ncvTest(lmauc)


## synID
auc_fold_change %>% 
    filter(type == "aggregate_fraction") %>% 
    group_by(dpi) %>% 
    wilcox_test(log2FC ~ synID, p.adjust.method = "holm")

## dpi
auc_fold_change %>% 
    filter(type == "aggregate_fraction" & log2FC != -Inf) %>% 
    wilcox_test(log2FC ~ dpi, detailed = TRUE)

## strain
auc_fold_change %>% 
    filter(type == "aggregate_fraction") %>% 
    group_by(dpi) %>% 
    kruskal_test(log2FC ~ strain)

auc_fold_change %>% 
    filter(type == "aggregate_fraction") %>% 
    group_by(dpi) %>% 
    dunn_test(log2FC ~ strain, p.adjust.method = "holm")
