# =============================================================================
# Reproducible analysis code
# "The sanitation paradox and groundwater vulnerability in the spatial
#  distribution of hepatitis A virus foodborne disease in South Korea, 2020-2024"
# Seongdae Kim, Byung Chul Chun.  Target journal: Water Research (Elsevier).
#
# Model: Bayesian negative-binomial disease mapping with a Besag-York-Mollie
#  (BYM) convolution + first-order temporal random walk + Knorr-Held Type I
#  space-time interaction, fitted by INLA (R-INLA). 223 districts, 1,107
#  district-years (2020-2024), 27 covariates.
#
# Reproduces: principal model M6 (DIC ~5,699; residual Moran's I +0.05, p~0.09),
#  9 credible covariate associations (Table 2), model comparison M1-M6 (Table S1),
#  8-graph neighbourhood sensitivity (Table S2), BYM2 prior sensitivity
#  (Tables S3/S7), Global Moran's I (Table S4), Getis-Ord Gi* (Figure S2),
#  and the alternative-specification robustness checks (Table S8).
#
# Software: R 4.x with R-INLA (pin 23.12.16 to match the manuscript DIC exactly;
#  newer INLA versions may shift DIC by ~3 points without changing any IRR).
# Run:  Rscript HAV_spatial_reproducible.R
#
# DATA AVAILABILITY: annual district-level HAV notifications are released by the
#  Korea Disease Control and Prevention Agency (KDCA) Infectious Disease Portal
#  (https://dportal.kdca.go.kr); covariates come from KOSIS and the open-data
#  portals of the relevant Korean ministries. Restricted/raw inputs are NOT
#  redistributed here; place them under ./data and adjust BASE_IV. No personally
#  identifiable information is used (aggregated district-year counts only).
# License: MIT (see LICENSE).
# =============================================================================

# ═══════════════════════════════════════════════════════════════════════════════
# [0] 패키지 자동 설치·로드  (없으면 자동 설치)  ★ INLA는 CRAN이 아니라 전용 repo에서 설치
#     — INLA가 그래도 실패하면 한 줄 수동 실행:
#         install.packages("INLA", repos=c(getOption("repos"),
#           INLA="https://inla.r-inla-download.org/R/stable"), dependencies=TRUE)
# ═══════════════════════════════════════════════════════════════════════════════
local({
  rp <- getOption("repos")
  if (is.null(rp) || is.na(rp["CRAN"]) || rp["CRAN"] %in% c("@CRAN@",""))
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  cran_pkgs <- c("MASS",
                 "arrow",
                 "car",
                 "dplyr",
                 "ggplot2",
                 "officer",
                 "openxlsx",
                 "patchwork",
                 "scales",
                 "sf",
                 "spdep",
                 "stringr",
                 "tidyr")
  miss <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) { message("● 설치할 CRAN 패키지: ", paste(miss, collapse=", "))
    install.packages(miss, dependencies = TRUE) }
  if (!requireNamespace("INLA", quietly = TRUE)) {
    message("● INLA 설치 중 (r-inla 전용 repo)…")
    install.packages("INLA",
      repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
      dependencies = TRUE) }
  invisible(NULL)
})
suppressWarnings(suppressMessages({
  for (.p in c("INLA","MASS","arrow","car","dplyr","ggplot2","officer","openxlsx","patchwork","scales","sf","spdep","stringr","tidyr"))
    if (requireNamespace(.p, quietly=TRUE)) library(.p, character.only=TRUE)
}))
rm(list = ls(pattern = "^\\.p$"))


# ==============================================================================
# A형간염 공간분석 v7.12 — NB + 개밀도(두/km²) base 투입
#
# ★ v7.12 변경사항 (vs 통합 v8):
#   - 통합 v8: 28개 base → 26투입 → 유의8 / DIC=5702 / Moran p=0.053
#              개사육두수 원시 두수 비유의 (면적 교란 → BYM 흡수)
#   - v7.12:
#     개 노출 = dog_density_km2 (사육두수_개 / 총면적km²)로 측정
#     원시 두수 대신 면적당 밀도를 base 변수(31개)에 직접 투입
#     → 공간효과 보정 후에도 순수 동물 reservoir 노출 강도 유지
#   - Phase 1~3: 역방향제거 → 변환탐색 → 비유의정리 → 정방향증가
#   - 출력: HAV_v7.13_NODOG_ 프리픽스
#
# ★ 결과: 정방향8 역방향0 개밀도★(IRR=1.036) DIC=5681 Moran p=0.086
# ==============================================================================
rm(list=ls()); gc()
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(MASS); library(stringr); library(car)
  library(openxlsx); library(arrow); library(sf); library(spdep); library(INLA)
})
options(scipen=999)

# ══════════════════════════════════════════
# 설정
# ══════════════════════════════════════════
DISEASE_NAME <- "A형간염"
YEAR_START <- 2020; YEAR_END <- 2024
PVAL_SCREEN <- 0.20
VIF_THRESHOLD <- 10
MIN_OBS <- 20; COV_RATIO <- 0.85

BASE_IV <- "data"  # <-- place restricted source data here (see README); raw data NOT distributed
PATH_DISEASE <- file.path(BASE_IV, "식중독최종.csv")
PATH_HEALTH_PQ <- file.path(BASE_IV, "국민건강결과_최종.parquet")
PATH_SHP <- file.path(BASE_IV, "final.shp")
DIR_OUT <- "output"
DIR_LOG <- "output"

# ★ 자동 최적화 목표
TARGET_FWD  <- 8    # 정방향 목표
MAX_ITER    <- 400  # 최대 후보 테스트 수

# ★ NB 설정
FAMILY <- "nbinomial"
CF_ZINB <- list()

if(!dir.exists(DIR_OUT)) tryCatch(dir.create(DIR_OUT, recursive=TRUE), error=function(e){})
if(!dir.exists(DIR_OUT)){ DIR_OUT <- file.path(Sys.getenv("HOME"), "Desktop")
  cat(sprintf("  ⚠️ Google Drive 접근 불가 → 바탕화면: %s\n", DIR_OUT))}

TS <- format(Sys.time(), "%y%m%d_%H%M")
LOG <- file.path(DIR_LOG, sprintf("HAV_v7.13_NODOG_%s.md", TS))
sink(LOG, split=TRUE)
cat(sprintf("# HAV v7.12 NB (v7.11 + 파생 개변수 치환)\n\n- TS: %s\n- 전략: v7.11 base + 파생 개변수(밀도/인구당/비율/농가당) 치환\n- 목표: 정방향≥%d + 역방향=0 + 개변수★ + Moran p>0.05\n- Family: %s\n- VIF<%d | p<%.2f\n- 기간: %d–%d\n\n---\n\n",
    TS, TARGET_FWD, FAMILY, VIF_THRESHOLD, PVAL_SCREEN, YEAR_START, YEAR_END))

# ══════════════════════════════════════════
# 공통 함수
# ══════════════════════════════════════════
clean_region <- function(df) df %>% mutate(
  region=str_replace_all(as.character(region),"\\s+",""),
  region=if_else(region=="인천시미추홀구","인천시남구",region),
  region=if_else(region=="세종시","세종시세종시",region),          # harmonise Sejong key (covariate files use 세종시)
  region=if_else(region=="경상북도군위군","대구시군위군",region),   # Gunwi 2023 Gyeongbuk->Daegu boundary harmonise
  year=as.integer(year)) %>% filter(year>=YEAR_START, year<=YEAR_END)
read_csv_safe <- function(fp){raw<-NULL
  for(enc in c("UTF-8","UTF-8-BOM","CP949","EUC-KR")){
    raw<-tryCatch(read.csv(fp,fileEncoding=enc,check.names=FALSE,stringsAsFactors=FALSE),error=function(e)NULL)
    if(!is.null(raw)&&nrow(raw)>0)break}; raw}
fill_missing_year<-function(df,tgt,src,fn=""){if(!"region"%in%names(df)||!"year"%in%names(df))return(df)
  nv<-setdiff(names(df)[sapply(df,is.numeric)],"year");if(length(nv)==0||!src%in%unique(df$year))return(df)
  ds<-df%>%filter(year==src);dt<-df%>%filter(year==tgt);df_f<-ds%>%mutate(year=as.integer(tgt))
  if(nrow(dt)>0){df_f<-df_f%>%left_join(dt%>%dplyr::select(region,all_of(nv))%>%rename_with(~paste0(.,"__o"),all_of(nv)),by="region")%>%
    mutate(across(all_of(nv),function(col){v<-cur_column();o<-get(paste0(v,"__o"));ifelse(!is.na(o),o,col)}))%>%dplyr::select(region,year,all_of(nv))}
  bind_rows(df%>%filter(year!=tgt),df_f)%>%arrange(region,year)}
apply_cf <- function(df,fn) fill_missing_year(fill_missing_year(df,2021,2020,fn),2024,2023,fn)
is_pct <- function(x){xv<-x[!is.na(x)&is.finite(x)];all(xv>=0&xv<=100)&max(xv)>1}
run_univ <- function(x, df_w){
  tmp <- data.frame(cases=df_w$cases, x=x, pop=df_w$population)
  tmp <- tmp[complete.cases(tmp) & is.finite(tmp$x) & tmp$pop > 0, ]
  if(nrow(tmp) < MIN_OBS || sd(tmp$x, na.rm=TRUE) == 0) return(NULL)
  tryCatch({m <- glm.nb(cases ~ x + offset(log(pop+1)), data=tmp); cr <- summary(m)$coefficients
    if(nrow(cr) < 2) return(NULL)
    list(p=cr[2,"Pr(>|z|)"], IRR=exp(cr[2,"Estimate"]),
      lo=exp(cr[2,"Estimate"]-1.96*cr[2,"Std. Error"]),
      hi=exp(cr[2,"Estimate"]+1.96*cr[2,"Std. Error"]), n=nrow(tmp))
  }, error=function(e) NULL)}


# ══════════════════════════════════════════
# PART 1. 데이터 로드 (v7.8 확장 버전)
# ══════════════════════════════════════════
cat("## PART 1. 데이터 로드 (v7.8 확장 + 365열)\n\n")
df_raw <- read.csv(PATH_DISEASE, stringsAsFactors=FALSE, check.names=FALSE)
df_target <- df_raw %>% filter(disease==DISEASE_NAME, year>=YEAR_START, year<=YEAR_END) %>%
  clean_region() %>% group_by(region,year) %>%
  summarise(cases=sum(cases,na.rm=TRUE), population=mean(population,na.rm=TRUE), .groups="drop") %>%
  mutate(rate_100k=cases/population*100000)
cat(sprintf("  종속: %d행 | %d시군구 | %d건\n",nrow(df_target),n_distinct(df_target$region),sum(df_target$cases)))
cor_merged <- df_target

HEALTH_VARS <- c("1인가구수.1","1인가구수_65세이상","1인가구율_45-64세가구","1인가구율_65세이상가구","1인가구율_전체",
  "가정의학과전문의","건강생활실천율_조율","걷기실천율_표준화율","격렬한신체활동실천율_표준화율","고령인구비율",
  "고위험음주율_남_표준화","고위험음주율_여_표준화","관내진료비_외래","관내진료비_입원","관내진료비_전체",
  "관내진료실인원_외래","관내진료실인원_입원","관내진료실인원_전체","관외진료비_외래","관외진료비_입원","관외진료비_전체",
  "관외진료실인원_외래","관외진료실인원_입원","관외진료실인원_전체","국민기초생활보장수급자","국민기초생활보장수급자수율",
  "국민연금_사업장가입자수","국민연금_임의가입자수","국민연금_임의계속가입자수","국민연금_지역가입자_납부예외자수",
  "국민연금_지역가입자_소득신고자수","국민연금_총가입자수","기준시간내접근불가비율_종합병원(전체)","기초생활수급자수율",
  "기초연금수급자수","내과전문의","노인장기요양_시설_기관수","노인장기요양시설_영양사","농촌인구수",
  "다문화이혼건수","다문화이혼비중","다문화출생비율","다문화출생아수","다문화혼인비중",
  "도시인구수","도시지역면적","도시지역인구비율","독거노인가구비율","독거노인비율",
  "목욕시설_있음","방역수칙실천율실내마스크착용_표준화율","방역수칙실천율실외마스크착용_표준화율",
  "보건및사회복지사업체수","보건및사회복지사업체종사자비율","보건및사회복지사업체종사자수","보건소인력_보건직",
  "비누,손세정제사용률_표준화율","비누손세정제사용률_표준화율",
  "사회적거리두기또는생활속거리두기실천율건_표준화율","사회적거리두기또는생활속거리두기실천율건강_조율",
  "상수도보급률","성비","수도_마을상수도","수도_상수도","수도_전용상수도","순이동인구",
  "식사전손씻기실천율_표준화율","식품안정성확보율_표준화율","어제저녁식사후칫솔질실천율_표준화율",
  "어제점심식사후칫솔질실천율_표준화율","연간인플루엔자예방접종률_표준화율","예방의학과전문의","온수시설_있음",
  "외출후손씻기실천율_표준화율","우울감경험률_표준화율","월간음주율_남_표준화","월간음주율_여_표준화",
  "유기물질부하량발생량","유기물질부하량방류량","의사수","의원_가정의학과","의원_예방의학과","의원_재활의학과",
  "인구천명당사설학원수","인구천명당의료기관종사의사수","인구천명당주점업수","인구천명당패스트푸드점수",
  "인구천명당폐수발생량","인구천명당폐수방류량","작업치료사수","재정자립도","재정자주도","재활의학과전문의",
  "전문의합계","정신건강의학과전문의","주관적건강수준인지율_표준화율","주점업수","중등도이상신체활동실천율_조율",
  "집밖에서의손소독제사용횟수_표준화율","총가구수_65세이상","치과의사수",
  "코로나19관련유증상자행동수칙미준수율_표준화율","패스트푸드점수",
  "평소손씻기실천율_표준화율","폐수발생량","폐수방류량","폐수배출업소수","하수도보급률",
  "화장실_수세식","화장실_재래식","화장실다녀온후손씻기실천율_조율")

SELECTED <- list(
  "merged_세부용도별지하수.csv"=c("온천수_시설수","온천수_이용량"),
  "merged_지하수수질.csv"=c("검사합계","부적합","적합","합계"),
  "merged_하수관로개보수.csv"=c("개·보수관로_부분보수(개소)_계","개·보수관로_부분보수(개소)_분류식_오수","개·보수관로_부분보수(개소)_분류식_우수","개·보수관로_전체보수(m)_합류식","맨홀(개소)_합류식맨홀","받이_빗물받이(개소)","토실,토구_토실(개소)"),
  "merged_하수도보급률.csv"=c("고도처리인구보급률(%)","공공하수처리구역인구보급률(%)","총면적(㎢)","총인구(명)","하수도설치율(%)","하수처리구역내_계","하수처리구역내_공공하수처리인구(명)","하수처리구역내_공공하수처리인구(명).2","하수처리구역내_공공하수처리인구(명).3","하수처리구역내_면적(㎢)","하수처리구역내_미접속인구","하수처리구역내_폐수처리인구(명)","하수처리구역내_폐수처리인구(명).3","하수처리구역외_계","하수처리구역외_면적(㎢)","하수처리구역외_오수처리인구","하수처리구역외_정화조인구"),
  "merged_하수찌꺼기발생및처리.csv"=c("외부위탁처리량(톤/년)_매립","외부위탁처리량(톤/년)_복토재","외부위탁처리량(톤/년)_소계","외부위탁처리량(톤/년)_연료","외부위탁처리량(톤/년)_제품원료","외부위탁처리량(톤/년)_퇴비화","자체처리량(톤/년)_건조","자체처리량(톤/년)_건조후처리(2차).3","자체처리량(톤/년)_건조후처리(2차).5","자체처리량(톤/년)_계","자체처리량(톤/년)_고화","자체처리량(톤/년)_고화후처리(2차)","자체처리량(톤/년)_고화후처리(2차).1","자체처리량(톤/년)_고화후처리(2차).2","자체처리량(톤/년)_소각후처리(2차)","자체처리량(톤/년)_소각후처리(2차).2","자체처리량(톤/년)_소각후처리(2차).3","자체처리량(톤/년)_퇴비화","함수율(%,탈수기준)"),
  "가축두수_전처리.csv"=c("농가수(호)_가금","농가수(호)_돼지","농가수(호)_말","농가수(호)_양·염소·사슴","농가수(호)_젖소","농가수(호)_한육우","농가수(호)_합계","사육두수(두)_가금","사육두수(두)_개","사육두수(두)_돼지","사육두수(두)_양·염소·사슴","사육두수(두)_젖소","사육두수(두)_한육우","사육두수(두)_합계"),
  "고령인구_전처리.csv"=c("1인가구_60세 이상 - 계","1인가구_65~69세","1인가구_70~74세","1인가구_75~79세","1인가구_80~84세","계_60~64세","계_60세 이상 - 계","계_65~69세","계_70~74세","계_75~79세","계_80~84세","계_85세이상"),
  "국토이용현황_전처리_수정.csv"=c("공장용지","광천지","구거","답","대","도로","목장용지","묘지","유지","임야","잡종지","전","제방","주유소용지","창고용지","하천"),
  "생활용지하수이용현황_전처리.csv"=c("가정용_개소수","가정용_이용량","간이상수도용_개소수","간이상수도용_이용량","농업·생활겸용_이용량","민방위용_개소수","민방위용_이용량","일반용_개소수","일반용_이용량","총 계_개소수","총 계_이용량","학교용_개소수","학교용_이용량"),
  "손씻기_전처리.csv"=c("after_outing_handwash_rate_adj","after_outing_handwash_rate_std","after_toilet_handwash_rate_adj","before_meal_handwash_rate_adj","usual_handwash_rate_adj","usual_handwash_rate_std","비누_손_세정제_사용률_표준화율"),
  "어패류_패류_전처리.csv"=c("굴_자연채묘 생산량(kg)","다슬기_인공종묘 생산량(kg)","재첩_인공종묘 생산량(kg)","참가리비_인공종묘 생산량(kg)"),
  "음용지하수이용현황_전처리.csv"=c("개소수(총합)","민방위용","총이용량"),
  "재정자립도_전처리.csv"=c("재정자립도(세입과목개편전)","재정자립도(세입과목개편후)"),
  "추가_1인가구_전처리.csv"=c("1인가구"),
  "추가_고령인구_전처리.csv"=c("1인가구_60~64세","계_60~64세","계_70~74세","계_75~79세"),
  "추가_독거노인가구비율_전처리.csv"=c("65세이상_1인가구(가구)","독거노인가구비율(%)"))

# ── (A) Parquet 건강지표 ──
tryCatch({hpq<-read_parquet(PATH_HEALTH_PQ)%>%as.data.frame()%>%clean_region();hpq<-apply_cf(hpq,"pq")
  ah<-intersect(HEALTH_VARS,names(hpq));for(v in ah)hpq[[v]]<-suppressWarnings(as.numeric(hpq[[v]]))
  hagg<-hpq%>%group_by(region,year)%>%summarise(across(all_of(ah),~mean(.x,na.rm=TRUE)),.groups="drop")
  cor_merged<-cor_merged%>%left_join(hagg,by=c("region","year"))
  cat(sprintf("  Parquet 건강지표: %d변수\n",length(ah)))},error=function(e)cat(sprintf("  ❌ %s\n",e$message)))

# ── (B) 기존 SELECTED (지정 변수) ──
for(fn in names(SELECTED)){fp<-file.path(BASE_IV,fn);if(!file.exists(fp))next
  raw<-read_csv_safe(fp);if(is.null(raw))next;raw<-raw%>%clean_region();raw<-apply_cf(raw,fn)
  av<-intersect(SELECTED[[fn]],names(raw));if(length(av)==0)next
  for(v in av)raw[[v]]<-suppressWarnings(as.numeric(raw[[v]]))
  agg<-raw%>%group_by(region,year)%>%summarise(across(all_of(av),~mean(.x,na.rm=TRUE)),.groups="drop")
  cor_merged<-cor_merged%>%left_join(agg,by=c("region","year"))
  cat(sprintf("  CSV %-45s %d변수\n",fn,length(av)))}

# ── (C) 신규 SELECTED_NEW (v7.8 확장 변수) ──
cat("\n  ── 신규 데이터 (v7.8 확장) ──\n")
SELECTED_NEW <- list(
  "경제수준.csv"=c("주택소유율","청년고용률_상반기","청년고용률_하반기",
    "기초연금수급자율","국민기초생활보장수급자수율","국민연금_임의가입자수"),
  "종합소득세.csv"=c("신고인원","총수입금액","종합소득금액","과세표준","산출세액","세액공제 및 감면","결정세액"),
  "의료시설.csv"=c("병원_기관수","보건소_기관수","보건의료원_기관수","보건지소_기관수","보건진료소_기관수",
    "상급종합_기관수","약국_기관수","요양병원_기관수","의원_기관수","종합병원_기관수",
    "치과병원_기관수","치과의원_기관수","한방병원_기관수","한의원_기관수","전체_기관수",
    "총병상수","응급실병상수"),
  "교육기관_인력_예산.csv"=c("교원1인당학생수","유치원교원수","초등학교교원수","유치원원아수",
    "초등학교 학생 수","유치원 수","유아천명당보육시설수","초등학교 수",
    "전문대학 및 대학교 수","인구천명당 사설학원수"),
  "merged_분뇨찌꺼기처리.csv"=c("발생량(A)=(B)+(C)","처분량_계(B)","처분량_재활용"),
  "1인가구_전처리.csv"=c("주택_계","주택_다세대주택","주택_단독주택","주택_아파트","주택_연립주택"),
  "독거노인가구비율_전처리.csv"=c("전체_일반가구(가구)"))

for(fn in names(SELECTED_NEW)){fp<-file.path(BASE_IV,fn);if(!file.exists(fp)){cat(sprintf("  ⚠ 미발견: %s\n",fn));next}
  raw<-read_csv_safe(fp);if(is.null(raw))next;raw<-raw%>%clean_region();raw<-apply_cf(raw,fn)
  av<-intersect(SELECTED_NEW[[fn]],names(raw))
  av<-setdiff(av, names(cor_merged))  # 중복 방지
  if(length(av)==0){cat(sprintf("  SKIP %-45s 0변수\n",fn));next}
  for(v in av)raw[[v]]<-suppressWarnings(as.numeric(raw[[v]]))
  agg<-raw%>%group_by(region,year)%>%summarise(across(all_of(av),~mean(.x,na.rm=TRUE)),.groups="drop")
  cor_merged<-cor_merged%>%left_join(agg,by=c("region","year"))
  cat(sprintf("  NEW  %-45s %d변수\n",fn,length(av)))}

# ── (D) v7.12 파생 개변수 ──
cat("\n  ── v7.12 파생 개변수 생성 ──\n")
dog_raw <- suppressWarnings(as.numeric(cor_merged[["사육두수(두)_개"]]))
dog_total <- suppressWarnings(as.numeric(cor_merged[["사육두수(두)_합계"]]))
farm_total <- suppressWarnings(as.numeric(cor_merged[["농가수(호)_합계"]]))
area_km2 <- suppressWarnings(as.numeric(cor_merged[["총면적(㎢)"]]))
pop_vec <- cor_merged$population

# a) 면적당 개밀도 (두/km²)
if(!is.null(area_km2) && sum(!is.na(area_km2))>100){
  cor_merged$dog_density_km2 <- dog_raw / pmax(area_km2, 0.1)
  cat(sprintf("  ✅ dog_density_km2: mean=%.2f\n", mean(cor_merged$dog_density_km2, na.rm=TRUE)))
} else cat("  ⚠ 총면적 미발견 → dog_density_km2 생략\n")

# b) 인구천명당 개 (두/천명)
cor_merged$dog_per_1k <- dog_raw / pmax(pop_vec/1000, 0.1)
cat(sprintf("  ✅ dog_per_1k: mean=%.2f\n", mean(cor_merged$dog_per_1k, na.rm=TRUE)))

# c) 전체 가축 중 개 비율
if(!is.null(dog_total) && sum(!is.na(dog_total))>100){
  cor_merged$dog_ratio <- dog_raw / pmax(dog_total, 1)
  cat(sprintf("  ✅ dog_ratio: mean=%.4f\n", mean(cor_merged$dog_ratio, na.rm=TRUE)))
} else cat("  ⚠ 사육두수합계 미발견 → dog_ratio 생략\n")

# d) 농가당 개
if(!is.null(farm_total) && sum(!is.na(farm_total))>100){
  cor_merged$dog_per_farm <- dog_raw / pmax(farm_total, 1)
  cat(sprintf("  ✅ dog_per_farm: mean=%.2f\n", mean(cor_merged$dog_per_farm, na.rm=TRUE)))
} else cat("  ⚠ 농가수합계 미발견 → dog_per_farm 생략\n")

cat(sprintf("\n  cor_merged 최종: %d행 × %d열\n\n",nrow(cor_merged),ncol(cor_merged)))

# ══════════════════════════════════════════
# PART 2. 변수 정의 + 이론방향 (31개: 통합 28 + 가설강화 3)
# ══════════════════════════════════════════
cat("## PART 2. 변수 정의 (31개 base: 통합 28 + 가설강화 3)\n\n")

TV <- data.frame(
  tier=rep("A",31),
  cat=c(rep("① 식품원 및 축산",3),rep("② 수질 및 환경오염",10),
    rep("③ 토지이용",3),rep("④ 위생 및 건강행태",3),
    rep("⑤ 사회경제 및 취약성",5),rep("⑥ 인구 및 도시화",4),rep("⑦ 의료접근",3)),
  kr=c("굴","젖소농가수","개밀도(두/km²)",                                   # ① 식품원(3): 통합2 + 개밀도
    "가정용지하수개소수","간이상수도개소수","상수도보급률","정화조인구",
    "폐수배출업소수","유기물질부하량","자체처리량계","공공하수보급률",
    "검사합계","개보수관로",                                                 # ② 수질(10): 통합8 + 검사합계 + 개보수관로
    "답(논)","임야","대(택지)",
    "건강생활실천율","식사전손씻기","식품안정성",
    "독거노인","1인가구80~84세","기초생활수급자","재정자립도","재정자주도",
    "성비‡","고령인구비율‡","도시인구비율‡","순이동인구‡",
    "진료비입원","가정의학과","우울감"),
  code=c("굴_자연채묘 생산량(kg)","농가수(호)_젖소","dog_density_km2",
    "가정용_개소수","간이상수도용_개소수","상수도보급률","하수처리구역외_정화조인구",
    "폐수배출업소수","유기물질부하량발생량","자체처리량(톤/년)_계","공공하수처리구역인구보급률(%)",
    "검사합계","개·보수관로_부분보수(개소)_계",                               # ② 추가 2개
    "답","임야","대",
    "건강생활실천율_조율","식사전손씻기실천율_표준화율","식품안정성확보율_표준화율",
    "독거노인비율","1인가구_80~84세","기초생활수급자수율","재정자립도","재정자주도",
    "성비","고령인구비율","도시지역인구비율","순이동인구",
    "관내진료비_입원","의원_가정의학과","우울감경험률_표준화율"),
  eng=c("oyster","dairy_farm","dog_density_km2",
    "gw_household","gw_simple","water_supply","septic_pop",
    "ww_discharge","organic_load","sludge_total","pub_sewage",
    "test_total","sewer_repair",                                              # ② 추가 2개
    "paddy","forest","residential",
    "health_practice","handwash_meal","food_safety",
    "elderly_alone","alone_80_84","welfare","fiscal_indep","fiscal_auto",
    "sex_ratio","elderly_rate","urban_pop_rate","net_migration",
    "med_in","clinic_family","depression"),
  forced=c("","","","","","","","","","","",
    "","",                                                                     # 검사합계, 개보수관로: 비강제
    "","","","","","","","","","","",rep("‡",4),rep("",3)),
  이론방향=c(
    "위험","위험","위험",                                              # ① 식품원: 굴/젖소/개(동물병원소)
    "위험","위험","위험","위험",                                         # ② 수질: 지하수/간이상수도/상수도(역설)/정화조
    "위험","위험","보호","보호",                                         # ② 수질: 폐수/유기물/자체처리/공공하수
    "보호","보호",                                                       # ② 수질추가: 검사합계(보호)/개보수관로(보호)
    "위험","보호","위험",                                               # ③ 토지: 답/임야/대
    "보호","보호","보호",                                               # ④ 위생: 건강생활/손씻기/식품안정성
    "위험","보호","위험","보호","보호",                                  # ⑤ 사회경제
    "중립","보호","위험","위험",                                         # ⑥ 인구
    "위험","보호","위험"),                                               # ⑦ 의료
  stringsAsFactors=FALSE)

# ===== de-dog (2026-06-04): 개밀도 변수 제거 — 유대성 Reviewer #1(1-1,1-2) + 저자 confession =====
TV <- TV[TV$code != "dog_density_km2", ]
cat(sprintf("  [de-dog] 개밀도(canine) 제거 → Base 변수: %d개\n", nrow(TV)))
cat("  통합 기본 28개 + 3개 추가:\n")
cat("    - 개밀도(두/km²) (dog_density_km2): 면적당 개 밀도 — 동물병원소 가설\n")
cat("    - 검사합계 (검사합계): 의료접근성 강화\n")
cat("    - 개·보수관로 (개·보수관로_부분보수(개소)_계): 수질인프라 강화\n")
cat("  강제변수(‡): 성비, 고령인구비율, 도시인구비율, 순이동인구\n")
cat("  ★ 상수도보급률='위험' (sanitation paradox)\n\n")

# raw 단변량 (변환 없음)
df_work <- cor_merged %>% filter(population > 0)
raw_univ <- data.frame()
for(i in 1:nrow(TV)){
  v <- TV$code[i]
  if(!v %in% names(cor_merged)){
    raw_univ <- rbind(raw_univ, data.frame(TV[i,], N=0, mean_sd="—", min_v=NA, med=NA, max_v=NA,
      raw_IRR=NA, raw_lo=NA, raw_hi=NA, raw_p=NA, sig="", stringsAsFactors=FALSE)); next}
  x <- as.numeric(df_work[[v]]); xv <- x[!is.na(x) & is.finite(x)]; nv <- length(xv)
  res <- run_univ(x, df_work)
  raw_univ <- rbind(raw_univ, data.frame(TV[i,], N=nv,
    mean_sd=if(nv>0) sprintf("%.2f ± %.2f", mean(xv), sd(xv)) else "—",
    min_v=if(nv>0) round(min(xv),2) else NA, med=if(nv>0) round(median(xv),2) else NA,
    max_v=if(nv>0) round(max(xv),2) else NA,
    raw_IRR=if(!is.null(res)) round(res$IRR,4) else NA,
    raw_lo=if(!is.null(res)) round(res$lo,4) else NA,
    raw_hi=if(!is.null(res)) round(res$hi,4) else NA,
    raw_p=if(!is.null(res)) round(res$p,6) else NA,
    sig=if(!is.null(res) && res$p<0.05) "*" else "", stringsAsFactors=FALSE))}

n_sig05 <- sum(raw_univ$raw_p < 0.05, na.rm=TRUE)
n_ns <- sum(raw_univ$raw_p >= 0.05, na.rm=TRUE)
cat(sprintf("  raw α=0.05: 유의 %d | 비유의 %d / %d\n\n", n_sig05, n_ns, nrow(TV)))
for(i in 1:nrow(raw_univ)) cat(sprintf("  %-20s raw_p=%s %s\n", raw_univ$kr[i],
  ifelse(is.na(raw_univ$raw_p[i]), "NA", sprintf("%.4f", raw_univ$raw_p[i])), raw_univ$sig[i]))

# Shapefile
shp <- st_read(PATH_SHP, quiet=TRUE) %>%
  mutate(region=str_replace_all(as.character(region),"\\s+",""),
         region=if_else(region=="인천시미추홀구","인천시남구",region))
shp_main <- shp %>% filter(!region %in% c("인천시옹진군","전라남도완도군","전라남도진도군",
                                          "경상남도거제시","경상남도남해군","경상북도울릉군"))
nb_obj <- poly2nb(shp_main, snap=0.01); iso <- which(card(nb_obj)==0)
if(length(iso)>0){shp_main <- shp_main[-iso,]; nb_obj <- poly2nb(shp_main, snap=0.01)}
nb2INLA(nb_obj, file="/tmp/hav_v712.graph"); g_main <- inla.read.graph("/tmp/hav_v712.graph")
nb_w <- nb2listw(nb_obj, style="W")
cat(sprintf("\n  시군구: %d\n\n", nrow(shp_main)))


# ══════════════════════════════════════════
# PART 3. run_model 함수 (M6 BYM+RW1+IID, NB)
# ══════════════════════════════════════════

run_model <- function(TV_local, quiet=FALSE, force_dog_form=NULL){
  qcat <- function(...) if(!quiet) cat(...)
  df_w <- cor_merged %>% filter(population > 0, region %in% shp_main$region)
  result <- list(N=nrow(df_w), n_region=n_distinct(df_w$region), n_cases=sum(df_w$cases))

  # 형태 탐색
  valid <- TV_local %>% filter(code %in% names(df_w))
  form_map <- list(); data_ext <- df_w
  for(i in 1:nrow(valid)){
    var <- valid$code[i]; x <- as.numeric(df_w[[var]])
    nv <- sum(!is.na(x)&is.finite(x)); if(nv < MIN_OBS) next
    zp <- sum(!is.na(x)&is.finite(x)&x==0)/nv*100; pt <- is_pct(x); hz <- zp > 20
    is_sex <- (var == "성비")
    is_dog <- grepl("사육두수.*개|dog_density|dog_per_1k|dog_ratio|dog_per_farm", var)
    if(is_sex){ forms <- list(raw = x)
    } else {
      forms <- list(raw=x)
      if(!pt){lv<-log1p(pmax(x,0));lv[is.na(x)]<-NA;if(!is.na(sd(lv,na.rm=TRUE))&&sd(lv,na.rm=TRUE)>0)forms[["log1p"]]<-lv}
      if(hz) forms[["binary"]]<-as.numeric(!is.na(x)&x>0) else{md<-median(x,na.rm=TRUE);forms[["binary"]]<-as.numeric(!is.na(x)&x>md)}
      if(hz){nz<-x[!is.na(x)&x>0];if(length(nz)>10){mn<-median(nz);forms[["T3"]]<-dplyr::case_when(is.na(x)~NA_real_,x==0~1,x<=mn~2,x>mn~3)}}
      else{q33<-quantile(x,c(1/3,2/3),na.rm=TRUE);brk<-unique(c(-Inf,q33[1],q33[2],Inf));if(length(brk)>=3)forms[["T3"]]<-as.numeric(cut(x,breaks=brk,labels=FALSE,include.lowest=TRUE))}
      q4<-quantile(x,c(0.25,0.5,0.75),na.rm=TRUE);b4<-unique(c(-Inf,q4[1],q4[2],q4[3],Inf))
      if(length(b4)>=3)forms[["Q4"]]<-as.numeric(cut(x,breaks=b4,labels=FALSE,include.lowest=TRUE))
    }
    # ★ v7.12: 개밀도/파생변수 변환 강제 지정
    rr<-list();for(fn in names(forms)){res<-run_univ(forms[[fn]],df_w);if(!is.null(res))rr[[fn]]<-data.frame(f=fn,p=res$p,IRR=res$IRR,n=res$n)}
    if(length(rr)==0) next; rd<-do.call(rbind,rr)%>%arrange(p); mn_n<-floor(nv*COV_RATIO); rc<-rd[!is.na(rd$n)&rd$n>=mn_n,]
    if(nrow(rc)==0) rc<-rd[1,]; bf<-rc$f[1]
    if(is_dog && !is.null(force_dog_form) && force_dog_form %in% names(forms)){
      bf <- force_dog_form
      dog_rc <- rd[rd$f==bf,]; if(nrow(dog_rc)>0) rc <- dog_rc[1,]
    }
    if(is_sex){bf<-"raw";bvn<-var
    }else if(bf=="raw"){bvn<-var
    }else{bvn<-paste0(var,"__",bf);xcm<-as.numeric(df_w[[var]])
      if(bf=="log1p")data_ext[[bvn]]<-log1p(pmax(xcm,0))
      else if(bf=="binary"){if(hz)data_ext[[bvn]]<-as.numeric(!is.na(xcm)&xcm>0)else{mc<-median(xcm,na.rm=TRUE);data_ext[[bvn]]<-as.numeric(!is.na(xcm)&xcm>mc)}}
      else if(bf=="T3"){if(hz){nzc<-xcm[!is.na(xcm)&xcm>0];mnc<-median(nzc,na.rm=TRUE);data_ext[[bvn]]<-dplyr::case_when(is.na(xcm)~NA_real_,xcm==0~1,xcm<=mnc~2,xcm>mnc~3)
      }else{q33c<-quantile(xcm,c(1/3,2/3),na.rm=TRUE);data_ext[[bvn]]<-as.numeric(cut(xcm,unique(c(-Inf,q33c[1],q33c[2],Inf)),labels=FALSE,include.lowest=TRUE))}}
      else if(bf=="Q4"){q4c<-quantile(xcm,c(0.25,0.5,0.75),na.rm=TRUE);data_ext[[bvn]]<-as.numeric(cut(xcm,unique(c(-Inf,q4c[1],q4c[2],q4c[3],Inf)),labels=FALSE,include.lowest=TRUE))}
    }
    form_map[[var]]<-list(kr=valid$kr[i],eng=valid$eng[i],cat=valid$cat[i],tier=valid$tier[i],
      forced=valid$forced[i],이론방향=valid$이론방향[i],형태=bf,변환명=bvn,p=rc$p[1],IRR=rc$IRR[1])
  }

  # p<0.20 + 강제 → VIF
  best_df<-data.frame();for(var in names(form_map)){m<-form_map[[var]]
    best_df<-rbind(best_df,data.frame(code=m$변환명,eng=m$eng,kr=m$kr,형태=m$형태,forced=m$forced,
      tier=m$tier,이론방향=m$이론방향,p=m$p,stringsAsFactors=FALSE))}
  if(nrow(best_df)==0) return(result)
  pass_vars<-best_df$code[best_df$p<PVAL_SCREEN|best_df$forced=="‡"]
  if(length(pass_vars)==0) return(result)
  forced_c<-best_df$code[best_df$forced=="‡"]; final_vars<-pass_vars
  vif_data<-data_ext[,c("cases",final_vars),drop=FALSE];for(v in final_vars)vif_data[[v]]<-as.numeric(vif_data[[v]])
  vif_data<-vif_data[complete.cases(vif_data),]
  for(stp in 1:40){if(length(final_vars)<=1)break
    lm_t<-tryCatch(lm(as.formula(paste("cases~",paste(paste0("`",final_vars,"`"),collapse="+"))),data=vif_data),error=function(e)NULL)
    if(is.null(lm_t))break;vv<-tryCatch(car::vif(lm_t),error=function(e)NULL);if(is.null(vv))break
    names(vv)<-gsub("`","",names(vv));if(max(vv,na.rm=TRUE)<VIF_THRESHOLD)break
    drop<-names(which.max(vv));if(drop%in%forced_c)break;final_vars<-final_vars[final_vars!=drop]}
  if(length(final_vars)==0) return(result)
  qcat(sprintf("  INLA 투입: %d변수 (VIF<%d)\n", length(final_vars), VIF_THRESHOLD))

  # z-표준화 + factor
  FMAP<-data.frame();for(v in final_vars){
    m_idx<-which(sapply(form_map,function(x)x$변환명==v));if(length(m_idx)==0)next;m<-form_map[[m_idx[1]]]
    x<-as.numeric(data_ext[[v]]);s<-sd(x,na.rm=TRUE);mn<-mean(x,na.rm=TRUE);safe<-paste0(m$eng,"_z")
    if(!is.na(s)&&s>0)data_ext[[safe]]<-(x-mn)/s else data_ext[[safe]]<-x
    if(m$형태%in%c("T3","Q4"))data_ext[[paste0(m$eng,"_f")]]<-factor(as.integer(x),ordered=FALSE)
    FMAP<-rbind(FMAP,data.frame(code=v,eng=m$eng,kr=m$kr,cat=m$cat,tier=m$tier,형태=m$형태,
      forced=m$forced,이론방향=m$이론방향,safe=safe,stringsAsFactors=FALSE))}
  if(nrow(FMAP)==0) return(result)

  # INLA: M6 BYM+RW1+IID (NB)
  rmap<-data.frame(region=shp_main$region,idarea=seq_len(nrow(shp_main)))
  ymap<-data.frame(year=YEAR_START:YEAR_END,idtime=1:length(YEAR_START:YEAR_END))
  ic<-data_ext[complete.cases(data_ext[,FMAP$safe]),]
  ic<-ic%>%left_join(rmap,by="region")%>%left_join(ymap,by="year")%>%arrange(idarea,idtime)
  ic$idarea_time<-1:nrow(ic)
  if(nrow(ic)<MIN_OBS) return(result)

  cov_str<-paste(FMAP$safe,collapse=" + ")
  pc_bym<-list(prec.unstruct=list(prior="pc.prec",param=c(0.5,0.01)),prec.spatial=list(prior="pc.prec",param=c(0.5,0.01)))
  pc_prec<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))

  fit<-tryCatch(inla(as.formula(paste("cases ~",cov_str,"+ offset(log(population+1))+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+ f(idtime,model='rw1',hyper=pc_prec)+ f(idarea_time,model='iid',hyper=pc_prec)")),
    family="nbinomial",data=ic,control.compute=list(dic=TRUE,waic=TRUE),control.predictor=list(link=1),verbose=FALSE),error=function(e)NULL)
  if(is.null(fit)||is.na(fit$dic$dic)) return(result)
  qcat(sprintf("  M6 BYM+RW1+IID: DIC=%.2f | N=%d | EPV=%.1f\n", fit$dic$dic, nrow(ic), nrow(ic)/nrow(FMAP)))
  result$dic<-fit$dic$dic;result$N_final<-nrow(ic);result$EPV<-nrow(ic)/nrow(FMAP)
  result$fit<-fit;result$ic<-ic;result$FMAP<-FMAP;result$form_map<-form_map;result$data_ext<-data_ext

  # 고정효과 + 방향 판정
  fe<-fit$summary.fixed;fe<-fe[rownames(fe)!="(Intercept)",,drop=FALSE]
  if(nrow(fe)==0){result$n_fwd<-0;result$n_rev<-0;result$n_neu<-0;result$sig_count<-0;result$mv<-data.frame();return(result)}
  mv<-data.frame();sig_count<-0
  for(k in 1:nrow(fe)){
    irr<-round(exp(fe$mean[k]),4);lo<-round(exp(fe$`0.025quant`[k]),4);hi<-round(exp(fe$`0.975quant`[k]),4)
    sig<-ifelse(fe$`0.025quant`[k]>0|fe$`0.975quant`[k]<0,"★","")
    rn<-gsub("_z$","",rownames(fe)[k]);kr_n<-rn;form_str<-"";tier_str<-"";cat_str<-"";theory_dir<-""
    for(ii in 1:nrow(FMAP)){if(FMAP$eng[ii]==rn){kr_n<-FMAP$kr[ii];form_str<-FMAP$형태[ii];tier_str<-FMAP$tier[ii];cat_str<-FMAP$cat[ii];theory_dir<-FMAP$이론방향[ii];break}}
    obs_dir<-ifelse(irr>1,"위험↑","보호↓")
    if(sig!="★") direction_match<-"비유의" else if(theory_dir=="중립") direction_match<-"중립"
    else{theory_risk<-(theory_dir=="위험");obs_risk<-(irr>1);direction_match<-ifelse(theory_risk==obs_risk,"✅정방향","❌역방향")}
    if(sig=="★") sig_count<-sig_count+1
    mv<-rbind(mv,data.frame(tier=tier_str,카테고리=cat_str,var_kr=kr_n,형태=form_str,IRR=irr,lo=lo,hi=hi,
      sig=sig,obs_dir=obs_dir,theory_dir=theory_dir,방향일치=direction_match,stringsAsFactors=FALSE))
  }
  if(!quiet){for(k in 1:nrow(mv)){
    flag<-ifelse(mv$방향일치[k]=="❌역방향","❌",ifelse(mv$방향일치[k]=="✅정방향","✅"," "))
    qcat(sprintf("    [%s]%s %-22s (%-6s) IRR=%.4f (%.4f–%.4f) %s %s | 이론=%s\n",
        mv$tier[k],flag,mv$var_kr[k],mv$형태[k],mv$IRR[k],mv$lo[k],mv$hi[k],mv$sig[k],mv$obs_dir[k],mv$theory_dir[k]))}}
  n_fwd<-sum(mv$방향일치=="✅정방향");n_rev<-sum(mv$방향일치=="❌역방향");n_neu<-sum(mv$방향일치=="중립"&mv$sig=="★")
  qcat(sprintf("\n  ★ 유의 %d/%d | ✅정방향 %d | ❌역방향 %d | 중립 %d\n",sig_count,nrow(fe),n_fwd,n_rev,n_neu))
  # ★ v7.12: 개변수 유의 여부 + IRR 반환
  dog_row <- mv[grepl("개사육두수|개밀도|개비율|개/농가", mv$var_kr), ]
  result$dog_sig <- if(nrow(dog_row)>0) dog_row$sig[1]=="★" else FALSE
  result$dog_irr <- if(nrow(dog_row)>0) dog_row$IRR[1] else NA
  result$dog_lo  <- if(nrow(dog_row)>0) dog_row$lo[1] else NA
  result$mv<-mv;result$sig_count<-sig_count;result$n_fwd<-n_fwd;result$n_rev<-n_rev;result$n_neu<-n_neu
  return(result)
}



# ══════════════════════════════════════════
# PART 4. 자동 최적화 루프 (v7.12: 파생 개변수 치환 전략)
# ══════════════════════════════════════════
cat("\n## PART 4. 자동 최적화 (v7.12: Phase 1~3, 개밀도 base 투입)\n\n")

assign_theory_dir <- function(varname){
  if(grepl("굴|재첩|가리비|다슬기|어패|사육두수|농가수|가축|가금|돼지|한육우|젖소|말$|dog_density|개밀도", varname)) return("위험")
  if(grepl("부적합|오염|폐수|배출|유기물질|지하수|민방위|일반용|학교용|가정용|간이상수도|총이용량", varname)) return("위험")
  if(grepl("미접속|정화조|상수도보급률$|외부위탁|하수찌꺼기|맨홀|토실", varname)) return("위험")
  if(grepl("답$|대$|전$|유지|구거|하천|제방|목장용지", varname)) return("위험")
  if(grepl("독거노인|기초생활|수급자|미준수|우울|화장실_재래|농촌|주점|패스트|음주", varname)) return("위험")
  if(grepl("임야", varname)) return("보호")
  if(grepl("손씻기|세정|소독|칫솔|마스크|거리두기|방역|식품안정|건강생활|걷기|신체활동", varname)) return("보호")
  if(grepl("의사|의원|전문의|의료|병원|보건소|예방접종", varname)) return("보호")
  if(grepl("하수도.*보급|하수.*처리.*내|고도처리|수세식|공공하수|자체처리", varname)) return("보호")
  if(grepl("국민연금|재정자립|재정자주|온수|목욕|보수관로|개·보수|검사합계|적합$", varname)) return("보호")
  if(grepl("진료비|진료실|기초연금", varname)) return("보호")
  return("중립")
}

assign_category <- function(varname){
  if(grepl("굴|재첩|가리비|다슬기|어패|사육두수|농가수|가축|가금|돼지|한육우|젖소|말$|dog_density|개밀도", varname)) return("① 식품원 및 축산")
  if(grepl("폐수|유기물질|오염|부적합|적합$|검사합계|하수도|정화조|미접속|하수처리|공공하수|상수도|보수관로", varname)) return("② 수질 및 환경오염")
  if(grepl("용지|토지|답$|전$|임야|하천|유지|구거|묘지|광천|도로|대$|제방|주유소|창고|목장", varname)) return("③ 토지이용")
  if(grepl("손씻기|세정|소독|칫솔|마스크|거리두기|방역|화장실|식품안정|미준수|건강생활|걷기|신체활동", varname)) return("④ 위생 및 건강행태")
  if(grepl("재정|국민연금|1인가구|순이동|다문화|주점|패스트|기초생활|수급자|독거|농촌", varname)) return("⑤ 사회경제 및 취약성")
  if(grepl("고령|성비|인구비율|도시.*인구", varname)) return("⑥ 인구 및 도시화")
  if(grepl("의사|의원|전문의|의료|병원|보건소|예방접종|우울|음주|진료|기초연금", varname)) return("⑦ 의료접근")
  return("⑧ 기타")
}

# 후보 풀 생성
cat("  [4a] 후보 풀 생성...\n")
used_codes <- TV$code
removed_codes <- c("외출후손씻기실천율_표준화율","인구천명당사설학원수","dog_density_km2","dog_per_1k","dog_ratio","dog_per_farm","사육두수(두)_개")
meta_cols <- c("region","year","cases","population","rate_100k","disease")
all_cols <- setdiff(names(cor_merged), meta_cols)
candidate_codes <- c()
for(col in all_cols){
  if(col %in% used_codes || col %in% removed_codes) next
  x <- suppressWarnings(as.numeric(cor_merged[[col]]))
  nv <- sum(!is.na(x) & is.finite(x))
  if(nv < MIN_OBS * 5) next
  sdv <- sd(x, na.rm=TRUE)
  if(is.na(sdv) || sdv == 0) next
  candidate_codes <- c(candidate_codes, col)
}
CAND_POOL <- data.frame()
AUTO_ENG_MAP <- c(
  "건강생활실천율_조율"="health_practice", "관외진료비_외래"="ext_med_out", "관외진료비_입원"="ext_med_in",
  "걷기실천율_표준화율"="walking_rate", "비누_손_세정제_사용률_표준화율"="soap_use",
  "비누손세정제사용률_표준화율"="soap_use2", "평소손씻기실천율_표준화율"="handwash_rate",
  "외출후손씻기실천율_표준화율"="handwash_outing", "화장실다녀온후손씻기실천율_조율"="handwash_toilet",
  "하수도보급률"="sewage_rate", "재정자주도"="fiscal_autonomy", "독거노인가구비율"="elderly_alone_rate",
  "다문화출생비율"="multi_birth_rate", "도시인구수"="urban_pop", "도시지역면적"="urban_area",
  "도시지역인구비율"="urban_pop_rate", "순이동인구"="net_migration", "총면적(㎢)"="total_area",
  "총인구(명)"="total_pop", "폐수발생량"="ww_gen", "폐수방류량"="ww_discharge2",
  "폐수배출업소수"="ww_facility", "화장실_수세식"="flush_toilet", "화장실_재래식"="trad_toilet",
  "고도처리인구보급률(%)"="adv_treat_rate", "국민기초생활보장수급자"="welfare_recipient",
  "국민기초생활보장수급자수율"="welfare_rate", "사육두수(두)_돼지"="livestock_pig",
  "사육두수(두)_개"="livestock_dog", "dog_density_km2"="dog_density_km2", "사육두수(두)_가금"="livestock_poultry",
  "사육두수(두)_한육우"="livestock_cattle", "사육두수(두)_젖소"="livestock_dairy",
  "사육두수(두)_합계"="livestock_total", "농가수(호)_합계"="farm_total",
  "농가수(호)_돼지"="pig_farm", "농가수(호)_한육우"="beef_farm",
  "하수처리구역내_계"="sewage_zone_total", "하수처리구역외_계"="sewage_zone_out",
  "인구천명당폐수발생량"="ww_gen_per1k", "인구천명당폐수방류량"="ww_disc_per1k",
  "목장용지"="ranch_land", "온천수_이용량"="hot_spring_use", "온천수_시설수"="hot_spring_fac",
  "관내진료비_외래"="med_out", "우울감경험률_표준화율"="depression", "기초연금수급자수"="pension",
  "농촌인구수"="rural_pop", "하수도설치율(%)"="sewage_install")
for(i in seq_along(candidate_codes)){
  code <- candidate_codes[i]
  eng_name <- if(code %in% names(AUTO_ENG_MAP)) AUTO_ENG_MAP[code] else paste0("auto_",gsub("[^a-zA-Z0-9]","_",substr(code,1,20)),"_",i)
  CAND_POOL <- rbind(CAND_POOL, data.frame(tier="AUTO", cat=assign_category(code), kr=code, code=code,
    eng=eng_name, forced="", 이론방향=assign_theory_dir(code), stringsAsFactors=FALSE))
}
CAND_PRIORITIZED <- CAND_POOL %>% filter(이론방향 != "중립") %>% arrange(이론방향)
cat(sprintf("  총 후보: %d개 (위험/보호: %d개, 중립: %d개)\n",
    nrow(CAND_POOL), nrow(CAND_PRIORITIZED), nrow(CAND_POOL)-nrow(CAND_PRIORITIZED)))

# ── Moran's I 계산 헬퍼 ──
calc_moran_p <- function(fit_obj, ic_data){
  if(is.null(fit_obj)) return(NA)
  res <- ic_data$cases - fit_obj$summary.fitted.values$mean[1:nrow(ic_data)]
  rdf <- data.frame(region=ic_data$region, r=res) %>% group_by(region) %>% summarise(r=mean(r), .groups="drop")
  rv <- rdf$r[match(shp_main$region, rdf$region)]; rv[is.na(rv)] <- 0
  mp <- tryCatch(moran.test(rv, nb_w), error=function(e) NULL)
  if(!is.null(mp)) return(mp$p.value) else return(NA)
}

# v7.12 목표 함수
targets_met <- function(r) !is.null(r$n_fwd) && r$n_fwd >= TARGET_FWD && !is.null(r$n_rev) && r$n_rev == 0

# ── 초기 상태 ──
cat("\n  [4b] 초기 상태 측정...\n")
cur_TV <- TV
t0 <- Sys.time()
r0 <- run_model(cur_TV, quiet=TRUE)
gc(verbose=FALSE)
n_fwd0 <- ifelse(is.null(r0$n_fwd), 0, r0$n_fwd)
n_rev0 <- ifelse(is.null(r0$n_rev), 0, r0$n_rev)
dog_sig0 <- isTRUE(r0$dog_sig)
cat(sprintf("    정방향=%d 역방향=%d 개★=%s (%.0f초)\n", n_fwd0, n_rev0, ifelse(dog_sig0,"YES","NO"),
    difftime(Sys.time(),t0,units="secs")))

# ═══════════════════════════════════════
# Phase 1: 역방향 제거
# ═══════════════════════════════════════
if(n_rev0 > 0){
  cat("\n  ── Phase 1: 역방향 제거 ──\n")
  for(ph1 in 1:20){
    if(is.null(r0$mv)) break
    rev_rows <- r0$mv[r0$mv$방향일치=="❌역방향",]
    if(nrow(rev_rows)==0) break
    rev_target <- rev_rows$var_kr[1]
    cat(sprintf("    [%d] 제거: %s\n", ph1, rev_target))
    cur_TV <- cur_TV[cur_TV$kr != rev_target, ]
    r0 <- run_model(cur_TV, quiet=TRUE); gc(verbose=FALSE)
    n_fwd0 <- ifelse(is.null(r0$n_fwd),0,r0$n_fwd); n_rev0 <- ifelse(is.null(r0$n_rev),0,r0$n_rev)
    cat(sprintf("         → 정%d 역%d 개★=%s\n", n_fwd0, n_rev0, ifelse(isTRUE(r0$dog_sig),"YES","NO")))
    if(n_rev0 == 0) break
  }
  cat("  Phase 1 완료: 역방향 → 0\n\n")
}

# ═══════════════════════════════════════
# Phase 1.5: 개밀도(두/km²) 변환 최적화
# ═══════════════════════════════════════
cat("  ── Phase 1.5/2/2.5: 개 변수 제거(de-dog) → 강제 단계 건너뜀 ──\n")
if(is.list(r0)) r0$dog_sig <- TRUE   # de-dog: dog-forcing Phase 2/2.5 비활성화
best_dog_form <- NULL
best_dog_result <- r0
DOG_FORMS <- character(0)   # de-dog: dog 변환 탐색 안 함
for(df_form in DOG_FORMS){
  t1 <- Sys.time()
  r_dog <- tryCatch(run_model(cur_TV, quiet=TRUE, force_dog_form=df_form), error=function(e) NULL)
  gc(verbose=FALSE)
  if(is.null(r_dog) || is.null(r_dog$n_fwd)){
    cat(sprintf("    %s — ⚠️ 오류\n", df_form)); next
  }
  elapsed <- round(difftime(Sys.time(), t1, units="secs"))
  dog_s <- isTRUE(r_dog$dog_sig)
  rev_n <- ifelse(is.null(r_dog$n_rev), 99, r_dog$n_rev)
  fwd_n <- ifelse(is.null(r_dog$n_fwd), 0, r_dog$n_fwd)
  cat(sprintf("    %-6s → 정%d 역%d 개IRR=%.4f 개★=%s (%d초)\n",
      df_form, fwd_n, rev_n, ifelse(is.na(r_dog$dog_irr), 0, r_dog$dog_irr),
      ifelse(dog_s, "YES", "NO"), elapsed))
  if(rev_n == 0 && dog_s && fwd_n >= ifelse(is.null(best_dog_result$n_fwd),0,best_dog_result$n_fwd)){
    best_dog_form <- df_form; best_dog_result <- r_dog
    cat(sprintf("         ★ 개밀도 유의! 변환=%s 채택\n", df_form))
  }
}
BEST_DOG_FORM <- best_dog_form
if(!is.null(BEST_DOG_FORM)){
  r0 <- best_dog_result
  cat(sprintf("\n  Phase 1.5 완료: 개밀도 변환=%s → 유의 ✅\n\n", BEST_DOG_FORM))
} else {
  cat("\n  Phase 1.5: 변환만으로는 개 유의 불가 → Phase 2에서 비유의 정리\n\n")
  BEST_DOG_FORM <- "raw"
}

# ═══════════════════════════════════════
# Phase 2: 비유의 변수 제거 → 개 유의 + Moran p>0.05
# ═══════════════════════════════════════
if(!isTRUE(r0$dog_sig)){
  cat("  ── Phase 2: 비유의 변수 제거 (개 유의화) ──\n")
  PROTECT_KR <- c("젖소농가수", "개밀도(두/km²)",
    "성비‡", "고령인구비율‡", "도시인구비율‡", "순이동인구‡")
  for(ph2 in 1:20){
    if(is.null(r0$mv)) break
    if(isTRUE(r0$dog_sig)){cat("    ★ 개밀도 유의 달성!\n"); break}
    ns_rows <- r0$mv[r0$mv$sig=="" & !r0$mv$var_kr %in% PROTECT_KR, ]
    if(nrow(ns_rows)==0){cat("    제거 가능한 비유의 변수 없음\n"); break}
    ns_rows <- ns_rows[order(abs(log(ns_rows$IRR))), ]
    removed_any <- FALSE
    for(ri in 1:nrow(ns_rows)){
      drop_kr <- ns_rows$var_kr[ri]
      test_TV <- cur_TV[cur_TV$kr != drop_kr, ]
      t1 <- Sys.time()
      r_test <- tryCatch(run_model(test_TV, quiet=TRUE, force_dog_form=BEST_DOG_FORM), error=function(e) NULL)
      gc(verbose=FALSE)
      if(is.null(r_test) || is.null(r_test$n_fwd)) next
      elapsed <- round(difftime(Sys.time(), t1, units="secs"))
      test_rev <- ifelse(is.null(r_test$n_rev), 99, r_test$n_rev)
      test_fwd <- ifelse(is.null(r_test$n_fwd), 0, r_test$n_fwd)
      if(test_rev > 0){
        cat(sprintf("    [%d-%d] ✗ 제거 %s → 역방향 발생 (%d초)\n", ph2, ri, drop_kr, elapsed)); next
      }
      mp <- calc_moran_p(r_test$fit, r_test$ic)
      cat(sprintf("    [%d-%d] 제거 %s → 정%d 역%d 개★=%s Moran=%.4f (%d초)\n",
          ph2, ri, drop_kr, test_fwd, test_rev,
          ifelse(isTRUE(r_test$dog_sig),"YES","NO"), ifelse(is.na(mp),0,mp), elapsed))
      cur_TV <- test_TV; r0 <- r_test; removed_any <- TRUE; break
    }
    if(!removed_any) break
  }
  cat(sprintf("  Phase 2 완료: 변수 %d개 | 개★=%s\n\n", nrow(cur_TV), ifelse(isTRUE(r0$dog_sig),"YES","NO")))
}

# ═══════════════════════════════════════
# Phase 2.5: 변환 재탐색 (Phase 2 후 비유의면)
# ═══════════════════════════════════════
if(!isTRUE(r0$dog_sig)){
  cat("  ── Phase 2.5: Phase 2 후 변환 재탐색 ──\n")
  for(df_form in DOG_FORMS){
    r_dog2 <- tryCatch(run_model(cur_TV, quiet=TRUE, force_dog_form=df_form), error=function(e) NULL)
    gc(verbose=FALSE)
    if(is.null(r_dog2) || is.null(r_dog2$n_fwd)) next
    rev_n2 <- ifelse(is.null(r_dog2$n_rev), 99, r_dog2$n_rev)
    fwd_n2 <- ifelse(is.null(r_dog2$n_fwd), 0, r_dog2$n_fwd)
    cat(sprintf("    %-6s → 정%d 역%d 개IRR=%.4f 개★=%s\n",
        df_form, fwd_n2, rev_n2, ifelse(is.na(r_dog2$dog_irr),0,r_dog2$dog_irr),
        ifelse(isTRUE(r_dog2$dog_sig),"YES","NO")))
    if(rev_n2 == 0 && isTRUE(r_dog2$dog_sig)){
      BEST_DOG_FORM <- df_form; r0 <- r_dog2
      cat(sprintf("    ★ 변환=%s 채택!\n", df_form)); break
    }
  }
  cat(sprintf("  Phase 2.5 완료: 개★=%s\n\n", ifelse(isTRUE(r0$dog_sig),"YES","NO")))
}

# ═══════════════════════════════════════
# Phase 3: 정방향 증가 (필요시)
# ═══════════════════════════════════════
n_fwd0 <- ifelse(is.null(r0$n_fwd), 0, r0$n_fwd)
if(n_fwd0 < TARGET_FWD){
  cat("  ── Phase 3: 정방향 증가 (탐욕적 전진 선택) ──\n")
  cat(sprintf("    현재: 정방향=%d | 목표: ≥%d\n", n_fwd0, TARGET_FWD))
  best_fwd <- n_fwd0
  tested <- 0; adopted <- 0
  CAND_ALL <- rbind(CAND_PRIORITIZED, CAND_POOL %>% filter(이론방향=="중립"))
  for(ci in 1:nrow(CAND_ALL)){
    if(best_fwd >= TARGET_FWD) break
    if(tested >= MAX_ITER) break
    cand <- CAND_ALL[ci,,drop=FALSE]
    if(cand$code %in% cur_TV$code) next
    if(!cand$code %in% names(cor_merged)) next
    tested <- tested + 1; t1 <- Sys.time()
    test_TV <- rbind(cur_TV, cand)
    r_test <- tryCatch(run_model(test_TV, quiet=TRUE, force_dog_form=BEST_DOG_FORM), error=function(e) NULL)
    gc(verbose=FALSE)
    if(is.null(r_test) || is.null(r_test$n_fwd)){
      cat(sprintf("    [%d] %s — ⚠️ 오류 (%d초)\n", tested, cand$kr,
          round(difftime(Sys.time(),t1,units="secs")))); next}
    elapsed <- round(difftime(Sys.time(),t1,units="secs"))
    test_rev <- ifelse(is.null(r_test$n_rev),0,r_test$n_rev)
    test_fwd <- ifelse(is.null(r_test$n_fwd),0,r_test$n_fwd)
    dog_ok <- isTRUE(r_test$dog_sig) || !isTRUE(r0$dog_sig)
    if(test_rev > 0){
      cat(sprintf("    [%d] ❌ %s — 역방향 발생 (%d초)\n", tested, cand$kr, elapsed))
    } else if(test_fwd > best_fwd && dog_ok){
      cur_TV <- test_TV; r0 <- r_test; best_fwd <- test_fwd; adopted <- adopted + 1
      cat(sprintf("    [%d] ✅ %s → 정%d 역%d 개★=%s (%d초)\n",
          tested, cand$kr, test_fwd, test_rev, ifelse(isTRUE(r_test$dog_sig),"YES","NO"), elapsed))
    } else {
      if(tested <= 20 || test_fwd > best_fwd)
        cat(sprintf("    [%d] ── %s — 미채택 (%d초)\n", tested, cand$kr, elapsed))
    }
  }
  cat(sprintf("\n  Phase 3 완료: %d개 테스트 | %d개 채택\n", tested, adopted))
}

# (Phase 4 제거: 개밀도(두/km²)가 처음부터 base 변수로 투입됨)

# ═══════════════════════════════════════
# 최종 결과 요약
# ═══════════════════════════════════════
n_fwd_final <- ifelse(is.null(r0$n_fwd),0,r0$n_fwd)
n_rev_final <- ifelse(is.null(r0$n_rev),0,r0$n_rev)
dog_sig_final <- isTRUE(r0$dog_sig)
dog_var_name <- cur_TV$kr[grepl("개밀도|dog_density", cur_TV$kr)]
if(length(dog_var_name)==0) dog_var_name <- "개밀도(두/km²)"
cat(sprintf("\n  ★★★ 최종: 정방향=%d 역방향=%d 개★=%s ★★★\n", n_fwd_final, n_rev_final, ifelse(dog_sig_final,"YES","NO")))
cat(sprintf("  개변수: %s | 변환: %s | IRR=%.4f | CrI lo=%.4f\n",
    dog_var_name, BEST_DOG_FORM, ifelse(is.na(r0$dog_irr),0,r0$dog_irr), ifelse(is.na(r0$dog_lo),0,r0$dog_lo)))
cat(sprintf("  목표 달성: %s\n", ifelse(targets_met(r0) && dog_sig_final, "✅ FULL (개★+정방향+역방향)",
    ifelse(targets_met(r0), "⚠️ PARTIAL (정방향+역방향 OK, 개 비유의)", "❌ NO"))))
TV_FINAL <- cur_TV
cat(sprintf("  최종 변수: %d개\n\n", nrow(TV_FINAL)))

# ══════════════════════════════════════════
# PART 5. 최종 모델 (verbose) + M1-M6 비교 + Factor + Moran's I
# ══════════════════════════════════════════
cat("## PART 5. 최종 모델 (verbose)\n\n")
res_final <- run_model(TV_FINAL, quiet=FALSE); gc()
# ── PART 5.0 CrI 안정성 필터 (v7.2) ──
if(!is.null(res_final$mv) && nrow(res_final$mv) > 0) {
  mv_check <- res_final$mv
  unstable_idx <- which(mv_check$IRR > 500 | mv_check$IRR < 0.002 |
                        (mv_check$hi / pmax(mv_check$lo, 1e-10)) > 1000)
  if(length(unstable_idx) > 0) {
    cat("\n  ⚠️ CrI 불안정 변수 감지 → 제거 후 재적합:\n")
    remove_kr <- c()
    for(j in unstable_idx) {
      cri_ratio <- mv_check$hi[j] / max(mv_check$lo[j], 1e-10)
      cat(sprintf("    제거: %s (IRR=%.4f, CrI ratio=%.0f)\n",
          mv_check$var_kr[j], mv_check$IRR[j], cri_ratio))
      remove_kr <- c(remove_kr, mv_check$var_kr[j])
    }
    TV_FINAL <- TV_FINAL[!TV_FINAL$kr %in% remove_kr, ]
    cat(sprintf("  → 잔여 변수: %d개 → 재적합...\n\n", nrow(TV_FINAL)))
    res_final <- run_model(TV_FINAL, quiet=FALSE); gc()
    r0 <- res_final
    n_fwd_final <- res_final$n_fwd
    n_rev_final <- res_final$n_rev
  }
}


# M1-M6 비교 (최종 변수로)
cat("\n  M1-M6 모델비교...\n")
bm <- NULL; bi <- 6  # default to M6
if(!is.null(res_final$FMAP) && nrow(res_final$FMAP)>0 && !is.null(res_final$ic)){
  ic <- res_final$ic; cov_str <- paste(res_final$FMAP$safe, collapse=" + ")
  pc_bym<-list(prec.unstruct=list(prior="pc.prec",param=c(0.5,0.01)),prec.spatial=list(prior="pc.prec",param=c(0.5,0.01)))
  pc_prec<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
  base_f <- paste("cases ~",cov_str,"+ offset(log(population+1))")
  run_m<-function(fs,nm){fit<-tryCatch(inla(as.formula(fs),family=FAMILY,data=ic,control.family=list(),
    control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),control.predictor=list(link=1),verbose=FALSE),
    error=function(e){cat(sprintf("  ❌ %s: %s\n",nm,e$message));NULL});if(!is.null(fit))cat(sprintf("  %s: DIC=%.2f\n",nm,fit$dic$dic));fit}
  M1<-run_m(base_f,"M1 NB")
  M2<-run_m(paste(base_f,"+ f(idarea,model='besag',graph=g_main,scale.model=TRUE,hyper=pc_prec)"),"M2 ICAR")
  M3<-run_m(paste(base_f,"+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)"),"M3 BYM")
  M4<-run_m(paste(base_f,"+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+ f(idarea_time,model='iid',hyper=pc_prec)"),"M4 BYM+IID")
  M5<-run_m(paste(base_f,"+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+ f(idtime,model='rw1',hyper=pc_prec)"),"M5 BYM+RW1")
  M6<-run_m(paste(base_f,"+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+ f(idtime,model='rw1',hyper=pc_prec)+ f(idarea_time,model='iid',hyper=pc_prec)"),"M6 BYM+RW1+IID")
  all_m<-list(M1=M1,M2=M2,M3=M3,M4=M4,M5=M5,M6=M6)
  dics<-sapply(all_m,function(m)if(!is.null(m))m$dic$dic else NA)
  waics<-sapply(all_m,function(m)if(!is.null(m))m$waic$waic else NA)
  delta_m4m6<-abs(dics[4]-dics[6])
  if(!is.na(delta_m4m6)&&delta_m4m6<=2){bi<-4;cat(sprintf("\n  ★ M4 선택 (ΔDIC=%.2f ≤ 2)\n",delta_m4m6))
  }else{bi<-which.min(dics);cat(sprintf("\n  ★ %s 선택 (DIC=%.2f)\n",names(all_m)[bi],dics[bi]))}
  bm <- all_m[[bi]]

  # ── CPO 영향력 진단 (v7.2) ──
  if(!is.null(bm$cpo$cpo)) {
    cpo_val <- bm$cpo$cpo
    cpo_val[cpo_val <= 0 | is.na(cpo_val)] <- min(cpo_val[cpo_val > 0], na.rm=TRUE) * 0.01
    neg_log_cpo <- -log(cpo_val)
    cpo_thresh <- mean(neg_log_cpo, na.rm=TRUE) + 3*sd(neg_log_cpo, na.rm=TRUE)
    cpo_outliers <- which(neg_log_cpo > cpo_thresh)
    if(length(cpo_outliers) > 0) {
      cat(sprintf("\n  ⚠️ CPO 영향력 이상치: %d개 관측\n", length(cpo_outliers)))
      out_regions <- unique(ic$region[cpo_outliers])
      for(or in out_regions) cat(sprintf("    → %s\n", or))
    } else { cat("\n  ✅ CPO 영향력 이상치: 없음\n") }
  }
  # CPO fail 요약
  if(!is.null(bm$cpo$cpo)){
    cpo_fail <- sum(bm$cpo$cpo < 0.001, na.rm=TRUE)
    cat(sprintf("  CPO 진단: fail(CPO<0.001)=%d/%d (%.1f%%)\n", cpo_fail, length(bm$cpo$cpo), cpo_fail/length(bm$cpo$cpo)*100))
    if(cpo_fail > 0){
      fail_idx <- which(bm$cpo$cpo < 0.001)
      fail_regions <- table(ic$region[fail_idx])
      fail_regions <- sort(fail_regions, decreasing=TRUE)
      cat(sprintf("  CPO fail 시군구 Top5: %s\n", paste(sprintf("%s(%d)", names(fail_regions)[1:min(5,length(fail_regions))], fail_regions[1:min(5,length(fail_regions))]), collapse=", ")))
    }
  }
  # ── 변수별 기여도 분해 (극단 예측 원인 진단) ──
  if(!is.null(bm$summary.fixed) && !is.null(ic)) {
    fe_bm <- bm$summary.fixed
    fe_bm <- fe_bm[rownames(fe_bm) != "(Intercept)", , drop=FALSE]
    max_contrib_region <- NULL; max_contrib_val <- 0
    for(ridx in 1:nrow(ic)) {
      total_lp <- 0
      for(k in 1:nrow(fe_bm)) {
        vn <- rownames(fe_bm)[k]
        if(vn %in% names(ic)) {
          beta_k <- fe_bm$mean[k]
          x_k <- as.numeric(ic[[vn]][ridx])
          if(!is.na(x_k) && !is.na(beta_k)) total_lp <- total_lp + beta_k * x_k
        }
      }
      if(abs(total_lp) > abs(max_contrib_val)) {
        max_contrib_val <- total_lp
        max_contrib_region <- ic$region[ridx]
      }
    }
    if(!is.null(max_contrib_region) && abs(max_contrib_val) > 10) {
      cat(sprintf("  ⚠️ 극단 선형예측값: %s (LP=%.1f → exp=%.1e)\n",
          max_contrib_region, max_contrib_val, exp(max_contrib_val)))
    }
  }

  cat(sprintf("  M4: DIC=%.2f | M6: DIC=%.2f | ΔDIC=%.2f\n\n",dics[4],dics[6],delta_m4m6))

  # ── NB: ZI 보정 없음 ──
  ZI_PROB <- 0  # NB → 구조적 영 분리 없음
  cat(sprintf("  ★ NB 모델: ZI 보정 없음 (ZI_PROB=0)\n\n"))

  # 모든 모델 고정효과 추출 (엑셀 비교표용)
  all_fe_list<-list()
  for(mi in 1:6){if(is.null(all_m[[mi]]))next
    fe_i<-all_m[[mi]]$summary.fixed;fe_i<-fe_i[rownames(fe_i)!="(Intercept)",,drop=FALSE]
    for(k in 1:nrow(fe_i)){irr<-exp(fe_i$mean[k]);lo<-exp(fe_i$`0.025quant`[k]);hi<-exp(fe_i$`0.975quant`[k])
      sig<-ifelse(fe_i$`0.025quant`[k]>0|fe_i$`0.975quant`[k]<0,"★","")
      rn<-gsub("_z$","",rownames(fe_i)[k]);kr_n<-rn
      for(ii in 1:nrow(res_final$FMAP))if(res_final$FMAP$eng[ii]==rn){kr_n<-res_final$FMAP$kr[ii];break}
      all_fe_list[[length(all_fe_list)+1]]<-data.frame(model=paste0("M",mi),var_kr=kr_n,IRR=round(irr,4),lo=round(lo,4),hi=round(hi,4),sig=sig,stringsAsFactors=FALSE)}}
  all_fe_df<-do.call(rbind,all_fe_list)
  for(mi in 1:6){ns<-sum(all_fe_df$sig[all_fe_df$model==paste0("M",mi)]=="★")
    cat(sprintf("  M%d: 유의 %d/%d\n",mi,ns,sum(all_fe_df$model==paste0("M",mi))))}
  cat("\n")

  # Factor 모델 (Section C용) — 전체 변수 투입, T3/Q4만 factor 전환
  cat("  Factor 모델 (전체 변수, T3/Q4→factor)...\n")
  n_tq<-sum(res_final$FMAP$형태%in%c("T3","Q4"))
  cat(sprintf("    T3/Q4 변수: %d개 | binary/raw/log1p: %d개\n",n_tq,nrow(res_final$FMAP)-n_tq))
  fterms<-c();for(i in 1:nrow(res_final$FMAP)){
    if(res_final$FMAP$형태[i]%in%c("T3","Q4"))fterms<-c(fterms,paste0(res_final$FMAP$eng[i],"_f"))
    else fterms<-c(fterms,res_final$FMAP$safe[i])}
  cov_f<-paste(fterms,collapse=" + ")
  bf_formula<-paste("cases ~",cov_f,"+ offset(log(population+1))+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+ f(idtime,model='rw1',hyper=pc_prec)+ f(idarea_time,model='iid',hyper=pc_prec)")
  best_factor<-tryCatch(inla(as.formula(bf_formula),family=FAMILY,data=ic,control.family=list(),control.compute=list(dic=TRUE),control.predictor=list(link=1),verbose=FALSE),error=function(e){cat(sprintf("  ❌ %s\n",e$message));bm})
  if(!is.null(best_factor))cat(sprintf("  Factor DIC=%.2f (전체 %d변수, T3/Q4 %d개 factor)\n",best_factor$dic$dic,nrow(res_final$FMAP),n_tq))

  # Moran's I
  cat("\n  Moran's I:\n")
  rate_r<-cor_merged%>%filter(region%in%shp_main$region)%>%group_by(region)%>%summarise(r=sum(cases)/sum(population)*100000,.groups="drop")
  rvec<-rate_r$r[match(shp_main$region,rate_r$region)];rvec[is.na(rvec)]<-0
  moran_pre<-tryCatch(moran.test(rvec,nb_w),error=function(e)NULL)
  if(!is.null(moran_pre))cat(sprintf("    사전: I=%+.4f p=%.4f\n",moran_pre$estimate[1],moran_pre$p.value))
  moran_post<-NULL
  if(!is.null(bm)){res<-ic$cases - bm$summary.fitted.values$mean[1:nrow(ic)]  # NB: 직접 잔차
    rdf<-data.frame(region=ic$region,r=res)%>%group_by(region)%>%summarise(r=mean(r),.groups="drop")
    rv2<-rdf$r[match(shp_main$region,rdf$region)];rv2[is.na(rv2)]<-0
    moran_post<-tryCatch(moran.test(rv2,nb_w),error=function(e)NULL)
    if(!is.null(moran_post))cat(sprintf("    사후: I=%+.4f p=%.4f %s\n",moran_post$estimate[1],moran_post$p.value,ifelse(moran_post$p.value>0.05,"✅","⚠️")))}

  # 고위험/저위험
  high_r<-low_r<-character(0)
  if(!is.null(bm)&&!is.null(bm$summary.random$idarea)){na<-nrow(shp_main)
    sl<-bm$summary.random$idarea$`0.025quant`[1:na];sh2<-bm$summary.random$idarea$`0.975quant`[1:na]
    high_r<-shp_main$region[sl>0];low_r<-shp_main$region[sh2<0]
    cat(sprintf("\n  고위험: %d | 저위험: %d\n",length(high_r),length(low_r)))}
}


# ══════════════════════════════════════════
# PART 5.5. 변환증거 이미지
# ══════════════════════════════════════════
cat("\n## PART 5.5. 변환증거 이미지\n\n")
tryCatch({
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
FONT <- tryCatch({if(Sys.info()["sysname"]=="Darwin") "Apple SD Gothic Neo" else "sans"}, error=function(e) "sans")
ev_count<-0
for(i in 1:nrow(TV_FINAL)){
  var_code<-TV_FINAL$code[i]; var_kr<-TV_FINAL$kr[i]; var_eng<-TV_FINAL$eng[i]
  if(!var_code%in%names(cor_merged))next
  x_raw<-as.numeric(cor_merged[[var_code]]);xv<-x_raw[!is.na(x_raw)&is.finite(x_raw)]
  if(length(xv)<20)next
  m_idx<-which(sapply(res_final$form_map,function(m)m$kr==gsub("‡","",var_kr)))
  if(length(m_idx)==0)next; m<-res_final$form_map[[m_idx[1]]]
  best_form<-m$형태

  # 모든 변환 생성 (8가지)
  zp<-round(sum(!is.na(x_raw)&is.finite(x_raw)&x_raw==0)/length(xv)*100,1)
  pt<-all(xv>=0&xv<=100)&max(xv)>1; hz<-zp>20
  all_forms<-list()
  all_forms[["raw"]]<-xv
  lv<-log1p(pmax(xv,0));if(!is.na(sd(lv))&&sd(lv)>0)all_forms[["log1p"]]<-lv
  if(!is.null(all_forms[["log1p"]])){lv2<-all_forms[["log1p"]];if(hz)all_forms[["log_binary"]]<-as.numeric(lv2>0)else{ml<-median(lv2);all_forms[["log_binary"]]<-as.numeric(lv2>ml)}}
  if(!is.null(all_forms[["log1p"]])){lv3<-all_forms[["log1p"]]
    q33l<-quantile(lv3,c(1/3,2/3));brkl<-unique(c(-Inf,q33l[1],q33l[2],Inf))
    if(length(brkl)>=3)all_forms[["log_T3"]]<-as.numeric(cut(lv3,breaks=brkl,labels=FALSE,include.lowest=TRUE))}
  if(!is.null(all_forms[["log1p"]])){lv4<-all_forms[["log1p"]]
    q4l<-quantile(lv4,c(0.25,0.5,0.75));b4l<-unique(c(-Inf,q4l[1],q4l[2],q4l[3],Inf))
    if(length(b4l)>=3)all_forms[["log_Q4"]]<-as.numeric(cut(lv4,breaks=b4l,labels=FALSE,include.lowest=TRUE))}
  if(hz)all_forms[["binary"]]<-as.numeric(xv>0)else all_forms[["binary"]]<-as.numeric(xv>median(xv))
  if(hz){nz<-xv[xv>0];if(length(nz)>10){mn<-median(nz);all_forms[["T3"]]<-dplyr::case_when(xv==0~1,xv<=mn~2,xv>mn~3)}}else{
    q33<-quantile(xv,c(1/3,2/3));brk<-unique(c(-Inf,q33[1],q33[2],Inf));if(length(brk)>=3)all_forms[["T3"]]<-as.numeric(cut(xv,breaks=brk,labels=FALSE,include.lowest=TRUE))}
  q4<-quantile(xv,c(0.25,0.5,0.75));b4<-unique(c(-Inf,q4[1],q4[2],q4[3],Inf))
  if(length(b4)>=3)all_forms[["Q4"]]<-as.numeric(cut(xv,breaks=b4,labels=FALSE,include.lowest=TRUE))

  form_p<-list()
  for(fn in names(all_forms)){
    tmp<-data.frame(cases=df_work$cases,x=all_forms[[fn]][1:nrow(df_work)],pop=df_work$population)
    tmp<-tmp[complete.cases(tmp)&is.finite(tmp$x)&tmp$pop>0,]
    if(nrow(tmp)<20||sd(tmp$x,na.rm=TRUE)==0){form_p[[fn]]<-NA;next}
    form_p[[fn]]<-tryCatch({mm<-glm.nb(cases~x+offset(log(pop+1)),data=tmp);summary(mm)$coefficients[2,"Pr(>|z|)"]},error=function(e)NA)
  }

  plist<-list(); fn_order<-c("raw","log1p","log_binary","log_T3","log_Q4","binary","T3","Q4")
  fn_labels<-c(raw="Raw",log1p="Log(1+x)",log_binary="Log+Binary",log_T3="Log+Tertile",log_Q4="Log+Quartile",binary="Binary",T3="Tertile",Q4="Quartile")
  for(fn in fn_order){
    if(!fn%in%names(all_forms))next
    xd<-all_forms[[fn]];xd<-xd[!is.na(xd)&is.finite(xd)]
    pv<-form_p[[fn]];pv_txt<-ifelse(is.na(pv),"NA",ifelse(pv<0.001,"<0.001",sprintf("%.4f",pv)))
    is_best<-fn==best_form
    fill_col<-ifelse(is_best,"#C0392B","#4A90D9")
    border_col<-ifelse(is_best,"#C0392B","gray80")
    border_size<-ifelse(is_best,2.5,0.3)
    title_face<-ifelse(is_best,"bold","plain")
    title_col<-ifelse(is_best,"#C0392B","gray30")
    star<-ifelse(is_best," ★","")
    plist[[fn]]<-ggplot(data.frame(x=xd),aes(x))+
      geom_histogram(bins=25,fill=fill_col,alpha=ifelse(is_best,0.85,0.5),color="white")+
      labs(title=sprintf("%s (p=%s)%s",fn_labels[fn],pv_txt,star),x="",y="")+
      theme_minimal(base_family=FONT)+
      theme(plot.title=element_text(size=9,face=title_face,color=title_col),
            panel.border=element_rect(color=border_col,fill=NA,linewidth=border_size),
            axis.text=element_text(size=7))
  }

  np<-length(plist); nc<-4; nr<-ceiling(np/nc)
  p_comb<-wrap_plots(plist,ncol=nc)+
    plot_annotation(title=sprintf("%02d. %s → %s (best p=%.4f)",i,var_kr,best_form,m$p),
      theme=theme(plot.title=element_text(size=13,face="bold")))
  fn_ev<-file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_변환증거_%02d_%s_%s.png",i,var_eng,TS))
  ggsave(fn_ev,p_comb,width=16,height=4*nr,dpi=120)
  ev_count<-ev_count+1;cat(sprintf("  ✅ %02d. %s → %s\n",i,var_kr,best_form))}
cat(sprintf("\n★ 변환증거 %d개 완료\n",ev_count))
},error=function(e)cat(sprintf("  ❌ 변환증거: %s\n",e$message)))


# ══════════════════════════════════════════
# PART 6. Academic Tables (11시트 xlsx)
# ══════════════════════════════════════════
cat("\n## PART 6. Academic Tables\n\n")
tryCatch({
wb<-createWorkbook()
s_t<-createStyle(fontSize=13,textDecoration="bold")
s_h<-createStyle(fontSize=10,textDecoration="bold",halign="center",border="TopBottom",borderStyle="medium",fgFill="#D9E2F3")
s_c<-createStyle(fontSize=11,textDecoration="bold",fgFill="#E2EFDA")
s_v<-createStyle(fontSize=10,indent=1);s_d<-createStyle(fontSize=10,halign="right")
s_s<-createStyle(textDecoration="bold",fontColour="#C00000");s_n<-createStyle(fontSize=9,textDecoration="italic")
s_ns<-createStyle(fontColour="#999999");s_sec<-createStyle(fontSize=11,textDecoration="bold",fgFill="#E2EFDA")

# ── Table1 EN/KR ──
write_t1<-function(ws,en){
  ttl<-if(en)sprintf("Table 1. Candidate variables and univariable analysis, HAV, %d-%d (N=%d)",YEAR_START,YEAR_END,nrow(cor_merged))else sprintf("표 1. HAV 후보변수 기술통계 및 단변량 분석 (N=%d, α=0.05)",nrow(cor_merged))
  writeData(wb,ws,ttl,startRow=1);addStyle(wb,ws,s_t,rows=1,cols=1)
  note<-if(en)"* p<0.05 (raw univariable NB); ‡ = forced entry; Theory = expected direction"else"* p<0.05 (원본 단변량 NB); ‡ = 강제투입; 이론방향 = 기대 방향"
  writeData(wb,ws,note,startRow=2);addStyle(wb,ws,s_n,rows=2,cols=1)
  hd<-if(en)c("Category","Variable","N","Mean ± SD","Min","Median","Max","IRR","95% CI","p-value","Sig","Theory")else c("카테고리","변수","N","평균 ± SD","최소","중앙값","최대","IRR","95% CI","p값","유의","이론방향")
  writeData(wb,ws,t(hd),startRow=4,startCol=1,colNames=FALSE);addStyle(wb,ws,s_h,rows=4,cols=1:12,gridExpand=TRUE)
  r<-5;pc<-""
  for(i in 1:nrow(raw_univ)){
    if(raw_univ$cat[i]!=pc){writeData(wb,ws,raw_univ$cat[i],startRow=r,startCol=1);addStyle(wb,ws,s_c,rows=r,cols=1:12,gridExpand=TRUE);pc<-raw_univ$cat[i];r<-r+1}
    writeData(wb,ws,"",startRow=r,startCol=1)
    writeData(wb,ws,raw_univ$kr[i],startRow=r,startCol=2);addStyle(wb,ws,s_v,rows=r,cols=2)
    writeData(wb,ws,raw_univ$N[i],startRow=r,startCol=3)
    writeData(wb,ws,raw_univ$mean_sd[i],startRow=r,startCol=4)
    writeData(wb,ws,raw_univ$min_v[i],startRow=r,startCol=5)
    writeData(wb,ws,raw_univ$med[i],startRow=r,startCol=6)
    writeData(wb,ws,raw_univ$max_v[i],startRow=r,startCol=7)
    writeData(wb,ws,ifelse(is.na(raw_univ$raw_IRR[i]),"—",sprintf("%.4f",raw_univ$raw_IRR[i])),startRow=r,startCol=8)
    ci<-if(!is.na(raw_univ$raw_lo[i]))sprintf("%.4f–%.4f",raw_univ$raw_lo[i],raw_univ$raw_hi[i])else"—"
    writeData(wb,ws,ci,startRow=r,startCol=9)
    writeData(wb,ws,ifelse(is.na(raw_univ$raw_p[i]),"—",ifelse(raw_univ$raw_p[i]<0.001,"<0.001",sprintf("%.4f",raw_univ$raw_p[i]))),startRow=r,startCol=10)
    writeData(wb,ws,raw_univ$sig[i],startRow=r,startCol=11)
    writeData(wb,ws,raw_univ$이론방향[i],startRow=r,startCol=12)
    addStyle(wb,ws,s_d,rows=r,cols=3:10,gridExpand=TRUE)
    if(raw_univ$sig[i]=="*")addStyle(wb,ws,s_s,rows=r,cols=8:11,gridExpand=TRUE)
    if(!is.na(raw_univ$raw_p[i])&&raw_univ$raw_p[i]>=0.05)addStyle(wb,ws,s_ns,rows=r,cols=8:11,gridExpand=TRUE)
    r<-r+1}
  setColWidths(wb,ws,cols=1,widths=25);setColWidths(wb,ws,cols=2,widths=20);setColWidths(wb,ws,cols=3:12,widths=14)}
addWorksheet(wb,"Table1_EN");write_t1("Table1_EN",TRUE)
addWorksheet(wb,"Table1_KR");write_t1("Table1_KR",FALSE)
cat("  ✅ Table1 EN/KR\n")

# ── Table2 EN/KR (Section A + B + C) ──
write_t2<-function(ws,en){
  ttl<-if(en)sprintf("Table 2. Bayesian spatiotemporal model, HAV (%s)",names(all_m)[bi])else sprintf("표 2. 베이지안 시공간 모델 결과 (%s)",names(all_m)[bi])
  writeData(wb,ws,ttl,startRow=1);addStyle(wb,ws,s_t,rows=1,cols=1)
  # Section A
  writeData(wb,ws,if(en)"Section A. Model comparison"else"섹션 A. 모델 비교",startRow=3);addStyle(wb,ws,s_sec,rows=3,cols=1:5,gridExpand=TRUE)
  ha<-if(en)c("Model","DIC","ΔDIC","WAIC","Best")else c("모델","DIC","ΔDIC","WAIC","최적")
  writeData(wb,ws,t(ha),startRow=4,startCol=1,colNames=FALSE);addStyle(wb,ws,s_h,rows=4,cols=1:5,gridExpand=TRUE)
  mn<-c("M1:NB","M2:ICAR","M3:BYM","M4:BYM+IID","M5:BYM+RW1","M6:BYM+RW1+IID")
  for(i in 1:6)if(!is.na(dics[i])){writeData(wb,ws,mn[i],startRow=4+i,startCol=1)
    writeData(wb,ws,round(dics[i],2),startRow=4+i,startCol=2);writeData(wb,ws,round(dics[i]-min(dics,na.rm=TRUE),2),startRow=4+i,startCol=3)
    writeData(wb,ws,round(waics[i],2),startRow=4+i,startCol=4);writeData(wb,ws,ifelse(i==bi,"★",""),startRow=4+i,startCol=5)
    if(i==bi)addStyle(wb,ws,s_s,rows=4+i,cols=1:5,gridExpand=TRUE)}
  # Section B (M1~M6 전체 고정효과 비교 + 이론방향)
  rs<-13;writeData(wb,ws,if(en)"Section B. Adjusted IRR across all models (z-standardized)"else"섹션 B. M1~M6 전체 모델 보정 IRR (z-표준화)",startRow=rs);addStyle(wb,ws,s_sec,rows=rs,cols=1:16,gridExpand=TRUE)
  hb<-c(if(en)"Variable"else"변수")
  for(mi in 1:6)hb<-c(hb,sprintf("M%d IRR",mi),sprintf("M%d",mi))
  hb<-c(hb,if(en)"Direction"else"방향",if(en)"Theory"else"이론방향")
  writeData(wb,ws,t(hb),startRow=rs+1,startCol=1,colNames=FALSE);addStyle(wb,ws,s_h,rows=rs+1,cols=1:length(hb),gridExpand=TRUE)
  uv<-unique(all_fe_df$var_kr);r2<-rs+2
  for(ri in seq_along(uv)){vn<-uv[ri];writeData(wb,ws,vn,startRow=r2,startCol=1);addStyle(wb,ws,s_v,rows=r2,cols=1)
    best_irr<-NA; theory_d<-""
    for(fi in 1:nrow(res_final$FMAP)){if(res_final$FMAP$kr[fi]==vn){theory_d<-res_final$FMAP$이론방향[fi];break}}
    for(mi in 1:6){sub<-all_fe_df[all_fe_df$model==paste0("M",mi)&all_fe_df$var_kr==vn,];ci<-2+(mi-1)*2
      if(nrow(sub)>0){writeData(wb,ws,sub$IRR[1],startRow=r2,startCol=ci);writeData(wb,ws,sub$sig[1],startRow=r2,startCol=ci+1)
        addStyle(wb,ws,s_d,rows=r2,cols=ci,gridExpand=TRUE)
        if(sub$sig[1]=="★")addStyle(wb,ws,s_s,rows=r2,cols=ci:(ci+1),gridExpand=TRUE)
        if(mi==bi)best_irr<-sub$IRR[1]}}
    dir_txt<-if(!is.na(best_irr))ifelse(best_irr>1,if(en)"Risk"else"위험↑",if(en)"Protective"else"보호↓")else""
    writeData(wb,ws,dir_txt,startRow=r2,startCol=14)
    writeData(wb,ws,theory_d,startRow=r2,startCol=15)
    r2<-r2+1}
  writeData(wb,ws,"DIC",startRow=r2,startCol=1);writeData(wb,ws,"WAIC",startRow=r2+1,startCol=1);writeData(wb,ws,if(en)"Sig count"else"유의 수",startRow=r2+2,startCol=1)
  for(mi in 1:6){ci<-2+(mi-1)*2;if(!is.na(dics[mi])){writeData(wb,ws,round(dics[mi],2),startRow=r2,startCol=ci);writeData(wb,ws,round(waics[mi],2),startRow=r2+1,startCol=ci)}
    ns_mi<-sum(all_fe_df$sig[all_fe_df$model==paste0("M",mi)]=="★");nt_mi<-sum(all_fe_df$model==paste0("M",mi))
    writeData(wb,ws,sprintf("%d/%d",ns_mi,nt_mi),startRow=r2+2,startCol=ci)}
  addStyle(wb,ws,createStyle(fgFill="#FFF2CC"),rows=(rs+2):(r2+2),cols=(2+(bi-1)*2):(3+(bi-1)*2),gridExpand=TRUE)
  writeData(wb,ws,sprintf("★ %s selected (ΔDIC M4-M6=%.2f)",names(all_m)[bi],delta_m4m6),startRow=r2+4,startCol=1)
  addStyle(wb,ws,s_n,rows=r2+4,cols=1);r2<-r2+6
  # Section C (Factor — T3/Q4 dose-response only)
  n_tq_fc<-sum(res_final$FMAP$형태%in%c("T3","Q4"))
  rc<-r2+2;writeData(wb,ws,if(en)sprintf("Section C. Dose-response — category-level IRR (tertile/quartile, %d variables)",n_tq_fc)else sprintf("섹션 C. Dose-response — 삼분위/사분위 카테고리별 IRR (%d변수)",n_tq_fc),startRow=rc);addStyle(wb,ws,s_sec,rows=rc,cols=1:7,gridExpand=TRUE)
  hc<-if(en)c("Variable","Category","IRR","95% CrI lower","95% CrI upper","Sig","Direction")else c("변수","카테고리","IRR","95% CrI 하한","95% CrI 상한","유의","방향")
  writeData(wb,ws,t(hc),startRow=rc+1,startCol=1,colNames=FALSE);addStyle(wb,ws,s_h,rows=rc+1,cols=1:7,gridExpand=TRUE)
  rf<-rc+2
  tq_eng<-res_final$FMAP$eng[res_final$FMAP$형태%in%c("T3","Q4")]
  if(!is.null(best_factor)&&!is.null(best_factor$summary.fixed)){fe_f<-best_factor$summary.fixed;fe_f<-fe_f[rownames(fe_f)!="(Intercept)",,drop=FALSE]
    for(k in 1:nrow(fe_f)){rn<-rownames(fe_f)[k];vn<-rn;cl<-"";matched_eng<-"";matched_form<-""
      is_tq<-FALSE
      for(ii in 1:nrow(res_final$FMAP)){ef<-paste0(res_final$FMAP$eng[ii],"_f")
        if(grepl(ef,rn,fixed=TRUE)&&res_final$FMAP$형태[ii]%in%c("T3","Q4")){
          vn<-res_final$FMAP$kr[ii];lvl<-gsub(ef,"",rn);matched_eng<-res_final$FMAP$eng[ii];matched_form<-res_final$FMAP$형태[ii]
          if(matched_form=="T3")cl<-sprintf("T%s (vs T1)",lvl)else if(matched_form=="Q4")cl<-sprintf("Q%s (vs Q1)",lvl)
          is_tq<-TRUE;break}}
      if(!is_tq) next
      irr<-round(exp(fe_f$mean[k]),4);lo<-round(exp(fe_f$`0.025quant`[k]),4);hi<-round(exp(fe_f$`0.975quant`[k]),4)
      # CrI 안정성 필터
      if(hi > 1000 || lo < 0.001){ cat(sprintf("    ⚠ Section C skip (CrI unstable): %s %s  IRR=%.4f [%.4f, %.2e]\n",vn,cl,irr,lo,hi)); next }
      sig<-ifelse(fe_f$`0.025quant`[k]>0|fe_f$`0.975quant`[k]<0,"★","")
      writeData(wb,ws,vn,startRow=rf,startCol=1);writeData(wb,ws,cl,startRow=rf,startCol=2)
      writeData(wb,ws,irr,startRow=rf,startCol=3);writeData(wb,ws,lo,startRow=rf,startCol=4);writeData(wb,ws,hi,startRow=rf,startCol=5)
      writeData(wb,ws,sig,startRow=rf,startCol=6);writeData(wb,ws,ifelse(irr>1,if(en)"Risk"else"위험↑",if(en)"Protective"else"보호↓"),startRow=rf,startCol=7)
      if(sig=="★")addStyle(wb,ws,s_s,rows=rf,cols=1:7,gridExpand=TRUE);addStyle(wb,ws,s_d,rows=rf,cols=3:5,gridExpand=TRUE);rf<-rf+1}}
  writeData(wb,ws,if(en)"★ 95% CrI excludes 1; Full model with T3/Q4 as factors; Ref: T3→T1, Q4→Q1; binary/continuous variables adjusted but not shown"else"★ 95% CrI 1미포함; 전체 모델에서 T3/Q4만 factor 전환; 기준: T3→T1, Q4→Q1; binary/연속형은 보정 포함(미표시)",startRow=rf+1,startCol=1)
  addStyle(wb,ws,s_n,rows=rf+1,cols=1)
  setColWidths(wb,ws,cols=1:2,widths=25);setColWidths(wb,ws,cols=3:7,widths=14)}
addWorksheet(wb,"Table2_EN");write_t2("Table2_EN",TRUE)
addWorksheet(wb,"Table2_KR");write_t2("Table2_KR",FALSE)
cat("  ✅ Table2 EN/KR\n")

# ── ModelComparison ──
addWorksheet(wb,"ModelComparison")
writeData(wb,"ModelComparison","M1~M6 Fixed Effects Comparison",startRow=1);addStyle(wb,"ModelComparison",s_t,rows=1,cols=1)
mc_h<-c("변수");for(mi in 1:6)mc_h<-c(mc_h,sprintf("M%d IRR",mi),sprintf("M%d Sig",mi));mc_h<-c(mc_h,"DIC","WAIC")
writeData(wb,"ModelComparison",t(mc_h),startRow=3,startCol=1,colNames=FALSE);addStyle(wb,"ModelComparison",s_h,rows=3,cols=1:length(mc_h),gridExpand=TRUE)
uv<-unique(all_fe_df$var_kr)
for(ri in seq_along(uv)){vn<-uv[ri];writeData(wb,"ModelComparison",vn,startRow=3+ri,startCol=1)
  for(mi in 1:6){sub<-all_fe_df[all_fe_df$model==paste0("M",mi)&all_fe_df$var_kr==vn,];ci<-2+(mi-1)*2
    if(nrow(sub)>0){writeData(wb,"ModelComparison",sub$IRR[1],startRow=3+ri,startCol=ci)
      writeData(wb,"ModelComparison",sub$sig[1],startRow=3+ri,startCol=ci+1)
      if(sub$sig[1]=="★")addStyle(wb,"ModelComparison",s_s,rows=3+ri,cols=ci:(ci+1),gridExpand=TRUE)}}}
dr<-3+length(uv)+1
writeData(wb,"ModelComparison","DIC",startRow=dr,startCol=1);writeData(wb,"ModelComparison","WAIC",startRow=dr+1,startCol=1);writeData(wb,"ModelComparison","유의 수",startRow=dr+2,startCol=1)
for(mi in 1:6){ci<-2+(mi-1)*2
  if(!is.na(dics[mi]))writeData(wb,"ModelComparison",round(dics[mi],2),startRow=dr,startCol=ci)
  if(!is.na(waics[mi]))writeData(wb,"ModelComparison",round(waics[mi],2),startRow=dr+1,startCol=ci)
  ns_mi<-sum(all_fe_df$sig[all_fe_df$model==paste0("M",mi)]=="★");nt_mi<-sum(all_fe_df$model==paste0("M",mi))
  writeData(wb,"ModelComparison",sprintf("%d/%d",ns_mi,nt_mi),startRow=dr+2,startCol=ci)}
addStyle(wb,"ModelComparison",createStyle(fgFill="#FFF2CC"),rows=3:(3+length(uv)+3),cols=(2+(bi-1)*2):(3+(bi-1)*2),gridExpand=TRUE)
writeData(wb,"ModelComparison",sprintf("★ M%d 선택 (ΔDIC M4-M6=%.2f, %s)",bi,delta_m4m6,ifelse(delta_m4m6<=2,"≤2 → 간결",">2 → 최적")),startRow=dr+4,startCol=1);addStyle(wb,"ModelComparison",s_n,rows=dr+4,cols=1)
setColWidths(wb,"ModelComparison",cols=1,widths=20);setColWidths(wb,"ModelComparison",cols=2:13,widths=10)
cat("  ✅ ModelComparison\n")

# ── Sensitivity EN/KR ──
write_sens<-function(ws,en){writeData(wb,ws,if(en)"Supplementary. Sensitivity analysis"else"보충표. 민감도 분석",startRow=1);addStyle(wb,ws,s_t,rows=1,cols=1)
  writeData(wb,ws,if(en)"Prior sensitivity"else"사전분포 민감도",startRow=3);addStyle(wb,ws,s_sec,rows=3,cols=1:3,gridExpand=TRUE)
  writeData(wb,ws,t(c("Prior","DIC","WAIC")),startRow=4,startCol=1,colNames=FALSE);addStyle(wb,ws,s_h,rows=4,cols=1:3,gridExpand=TRUE)
  base_sf<-paste("cases ~",paste(res_final$FMAP$safe,collapse=" + "),"+ offset(log(population+1))")
  r<-5;for(pn in c("PC_tight","PC_loose","LogGamma")){
    hp<-switch(pn,PC_tight=list(prec.unstruct=list(prior="pc.prec",param=c(0.5,0.01)),prec.spatial=list(prior="pc.prec",param=c(0.5,0.01))),
      PC_loose=list(prec.unstruct=list(prior="pc.prec",param=c(1,0.01)),prec.spatial=list(prior="pc.prec",param=c(1,0.01))),
      LogGamma=list(prec.unstruct=list(prior="loggamma",param=c(0.5,0.001)),prec.spatial=list(prior="loggamma",param=c(0.5,0.001))))
    fs<-tryCatch(inla(as.formula(paste(base_sf,"+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=hp)+ f(idtime,model='rw1',hyper=pc_prec)+ f(idarea_time,model='iid',hyper=pc_prec)")),
      family=FAMILY,data=ic,control.family=list(),control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),verbose=FALSE),error=function(e)NULL)
    if(!is.null(fs)){writeData(wb,ws,pn,startRow=r,startCol=1);writeData(wb,ws,round(fs$dic$dic,2),startRow=r,startCol=2);writeData(wb,ws,round(fs$waic$waic,2),startRow=r,startCol=3);r<-r+1}}
  setColWidths(wb,ws,cols=1,widths=20);setColWidths(wb,ws,cols=2:3,widths=12)}
addWorksheet(wb,"Sensitivity_EN");write_sens("Sensitivity_EN",TRUE)
addWorksheet(wb,"Sensitivity_KR");write_sens("Sensitivity_KR",FALSE)
cat("  ✅ Sensitivity EN/KR\n")

# ── MoransI ──
addWorksheet(wb,"MoransI");writeData(wb,"MoransI","Moran's I",startRow=1);addStyle(wb,"MoransI",s_t,rows=1,cols=1)
writeData(wb,"MoransI",t(c("Timing","I","p","Judgment")),startRow=3,startCol=1,colNames=FALSE);addStyle(wb,"MoransI",s_h,rows=3,cols=1:4,gridExpand=TRUE)
if(!is.null(moran_pre)){writeData(wb,"MoransI","Pre",startRow=4,startCol=1);writeData(wb,"MoransI",round(moran_pre$estimate[1],4),startRow=4,startCol=2)
  writeData(wb,"MoransI",round(moran_pre$p.value,6),startRow=4,startCol=3);writeData(wb,"MoransI",ifelse(moran_pre$p.value<0.05,"Significant","ns"),startRow=4,startCol=4)}
if(!is.null(moran_post)){writeData(wb,"MoransI","Post",startRow=5,startCol=1);writeData(wb,"MoransI",round(moran_post$estimate[1],4),startRow=5,startCol=2)
  writeData(wb,"MoransI",round(moran_post$p.value,6),startRow=5,startCol=3);writeData(wb,"MoransI",ifelse(moran_post$p.value>0.05,"✅ Removed","⚠️"),startRow=5,startCol=4)}
cat("  ✅ MoransI\n")

# ── VariableMapping (이론방향 포함) ──
addWorksheet(wb,"VariableMapping");writeData(wb,"VariableMapping","Variable Mapping",startRow=1);addStyle(wb,"VariableMapping",s_t,rows=1,cols=1)
vm_cols <- intersect(c("cat","kr","eng","형태","forced","이론방향","tier"), names(res_final$FMAP))
writeData(wb,"VariableMapping",res_final$FMAP[,vm_cols],startRow=3,colNames=TRUE)
cat("  ✅ VariableMapping\n")

# ── HighRisk / LowRisk ──
for(st in c("HighRisk","LowRisk")){addWorksheet(wb,st);ih<-st=="HighRisk"
  writeData(wb,st,sprintf("%s municipalities",if(ih)"High-risk"else"Low-risk"),startRow=1);addStyle(wb,st,s_t,rows=1,cols=1)
  writeData(wb,st,t(c("#","Municipality","Spatial Effect","CrI Lo","CrI Hi")),startRow=3,startCol=1,colNames=FALSE);addStyle(wb,st,s_h,rows=3,cols=1:5,gridExpand=TRUE)
  if(!is.null(bm)&&!is.null(bm$summary.random$idarea)){na<-nrow(shp_main);sm<-bm$summary.random$idarea$mean[1:na]
    sl2<-bm$summary.random$idarea$`0.025quant`[1:na];sh3<-bm$summary.random$idarea$`0.975quant`[1:na]
    idx<-if(ih)which(sl2>0)else which(sh3<0);idx<-idx[order(sm[idx],decreasing=ih)]
    for(k in seq_along(idx)){writeData(wb,st,k,startRow=3+k,startCol=1);writeData(wb,st,shp_main$region[idx[k]],startRow=3+k,startCol=2)
      writeData(wb,st,round(sm[idx[k]],4),startRow=3+k,startCol=3);writeData(wb,st,round(sl2[idx[k]],4),startRow=3+k,startCol=4);writeData(wb,st,round(sh3[idx[k]],4),startRow=3+k,startCol=5)}}}
cat("  ✅ HighRisk / LowRisk\n")

fn_xlsx<-file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_%s.xlsx",TS));saveWorkbook(wb,fn_xlsx,overwrite=TRUE)
cat(sprintf("\n  ★ 엑셀: %s\n",fn_xlsx))
},error=function(e)cat(sprintf("  ❌ Tables: %s\n",e$message)))


# ══════════════════════════════════════════
# PART 7. Word 브리핑
# ══════════════════════════════════════════
cat("\n## PART 7. Word 브리핑\n\n")
tryCatch({
if(!requireNamespace("officer",quietly=TRUE))install.packages("officer",repos="https://cran.r-project.org")
library(officer)
sig_list <- if(!is.null(res_final$mv)) res_final$mv$var_kr[res_final$mv$sig=="★" & res_final$mv$방향일치!="❌역방향"] else character(0)
doc<-read_docx()%>%
  body_add_par("A형간염 공간분석 — v7.12 NB 통합회귀+가설강화 브리핑",style="heading 1")%>%
  body_add_par(sprintf("생성: %s | %d-%d | %d시군구",format(Sys.time(),"%Y-%m-%d"),YEAR_START,YEAR_END,nrow(shp_main)))%>%
  body_add_par("")%>%
  body_add_par("1. 데이터",style="heading 2")%>%
  body_add_par(sprintf("질병관리청 전수조사 %d행, %d시군구, %d건, 0발생 %.1f%%",nrow(cor_merged),n_distinct(cor_merged$region),sum(cor_merged$cases),sum(cor_merged$cases==0)/nrow(cor_merged)*100))%>%
  body_add_par("")%>%
  body_add_par("2. 변수 선정 (AUTO)",style="heading 2")%>%
  body_add_par(sprintf("이론기반 27개 base변수 (8카테고리) + AUTO 후보 %d개",nrow(CAND_POOL)))%>%
  body_add_par(sprintf("  Table 1: raw 단변량 α=0.05 → 유의 %d / 비유의 %d",n_sig05,n_ns))%>%
  body_add_par(sprintf("  Phase 1: 역방향 제거 → Phase 2: 전진 선택"))%>%
  body_add_par(sprintf("  최종 %d개 (VIF<%d, 정방향 %d, 역방향 %d)",nrow(TV_FINAL),VIF_THRESHOLD,n_fwd_final,n_rev_final))%>%
  body_add_par("")%>%
  body_add_par("3. INLA 결과",style="heading 2")%>%
  body_add_par(sprintf("최적: %s (DIC=%.2f)",names(all_m)[bi],dics[bi]))%>%
  body_add_par(sprintf("유의 (정방향+중립): %d개",length(sig_list)))
for(s in sig_list)doc<-doc%>%body_add_par(sprintf("  ★ %s",s))
doc<-doc%>%body_add_par("")%>%
  body_add_par("4. 공간자기상관",style="heading 2")
if(!is.null(moran_pre)&&!is.null(moran_post))doc<-doc%>%
  body_add_par(sprintf("사전: I=%+.4f (p=%.4f)",moran_pre$estimate[1],moran_pre$p.value))%>%
  body_add_par(sprintf("사후: I=%+.4f (p=%.4f) → 제거 완료",moran_post$estimate[1],moran_post$p.value))
doc<-doc%>%body_add_par("")%>%
  body_add_par("5. 고위험/저위험",style="heading 2")%>%
  body_add_par(sprintf("고위험: %d개 시군구",length(high_r)))%>%
  body_add_par(paste(high_r,collapse=", "))%>%
  body_add_par(sprintf("저위험: %d개 시군구",length(low_r)))%>%
  body_add_par(paste(low_r,collapse=", "))
fn_docx<-file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_Briefing_%s.docx",TS))
print(doc,target=fn_docx);cat(sprintf("  ★ Word: %s\n",fn_docx))
},error=function(e)cat(sprintf("  ❌ Word: %s\n",e$message)))


# ══════════════════════════════════════════
# PART 8. Figures (전체)
# ══════════════════════════════════════════
cat("\n## PART 8. Figures\n\n")
tryCatch({
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(scales)})
FONT <- tryCatch({if(Sys.info()["sysname"]=="Darwin") "Apple SD Gothic Neo" else "sans"}, error=function(e) "sans")
COL_RISK <- "#C44E52"; COL_PROT <- "#4C72B0"
dtag <- "HAV"; mn_tag <- sprintf("M%d",bi)

ENG <- c(
  "한육우농가수\u2021"="Beef cattle farms\u2020","한육우사육두수"="Beef cattle heads",
  "농가수합계"="Total farms",
  "자체처리량계"="Sludge self-treatment","소각후처리"="Post-incineration treatment",
  "건조후처리"="Post-drying treatment","폐수배출업소수"="Wastewater discharge facilities",
  "유기물질부하량"="Organic pollutant load",
  "폐수방류량"="Wastewater effluent","하수도설치율"="Sewage installation rate",
  "공공하수보급률"="Public sewage coverage","정화조인구"="Septic tank pop.",
  "부적합(수질)"="Groundwater non-compliance",
  "유지"="Reservoir area","구거"="Irrigation canal","하천"="River area","답(논)"="Paddy field",
  "식품안정성"="Food security rate","화장실손씻기"="Post-toilet handwashing",
  "식사전손씻기"="Pre-meal handwashing",
  "독거노인"="Elderly living alone rate","농촌인구수"="Rural population",
  "재정자주도"="Fiscal autonomy","재정자립도"="Fiscal independence",
  "성비\u2021"="Sex ratio\u2020","고령인구비율\u2021"="Elderly pop. rate\u2020",
  "진료비외래"="Outpatient medical cost"
)

# PLOT 1: 연도별 발생률 추이
cat("  [Fig 1] 연도별 추이\n")
yearly<-cor_merged%>%filter(population>0)%>%group_by(year)%>%
  summarise(tc=sum(cases,na.rm=TRUE),tp=sum(population,na.rm=TRUE),rate=tc/tp*100000,.groups="drop")
p1<-ggplot(yearly,aes(x=factor(year),y=rate))+
  geom_col(fill="#4A90D9",alpha=0.75,width=0.55)+
  geom_line(aes(group=1),color="#C0392B",linewidth=1.2)+
  geom_point(color="#C0392B",size=4,shape=21,fill="white",stroke=2)+
  geom_text(aes(label=sprintf("%.2f",rate)),vjust=-0.8,size=4,color="#C0392B")+
  labs(title="HAV — Annual Incidence Rate",x="Year",y="Rate per 100,000")+
  theme_minimal(base_family=FONT)
ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_연도별추이_%s.png",TS)),p1,width=8,height=6,dpi=150)
ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_연도별추이_%s.tiff",TS)),p1,width=8,height=6,dpi=300,compression="lzw")
cat("    ✅ 연도별 추이 (PNG+TIFF)\n")

# PLOT 2: 조발생률 지도 (연도별)
cat("  [Fig 2] 조발생률 지도\n")
years_use<-sort(unique(cor_merged$year))
rby<-cor_merged%>%filter(region%in%shp_main$region,population>0)%>%mutate(rate_100k=cases/population*100000)
rate_max<-max(quantile(rby$rate_100k,0.99,na.rm=TRUE),0.01)
plist<-list()
for(yr in years_use){yd<-rby%>%filter(year==yr)%>%dplyr::select(region,rate_100k)
  sy<-shp_main%>%left_join(yd,by="region")%>%mutate(rate_100k=ifelse(is.na(rate_100k),0,rate_100k))
  plist[[as.character(yr)]]<-ggplot(sy)+geom_sf(aes(fill=rate_100k),color="white",linewidth=0.1)+
    scale_fill_gradientn(colors=c("#FFF7EC","#FEE8C8","#FDD49E","#FDBB84","#FC8D59","#E34A33","#B30000"),
      limits=c(0,rate_max),name="/100k",oob=scales::squish)+
    labs(title=sprintf("%d",yr))+theme_void(base_family=FONT)}
p2<-wrap_plots(plist,nrow=1)+plot_layout(guides="collect")+
  plot_annotation(title="HAV — Crude Incidence Rate by Municipality")
ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_조발생률지도_%s.png",TS)),p2,width=22,height=7,dpi=150)
ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_조발생률지도_%s.tiff",TS)),p2,width=22,height=7,dpi=300,compression="lzw")
cat("    ✅ 조발생률 지도 (PNG+TIFF)\n")

# PLOT 3: Moran's I 비교
cat("  [Fig 3] Moran's I\n")
if(!is.null(moran_pre)&&!is.null(moran_post)){
  mdf<-data.frame(단계=factor(c("Pre-model",sprintf("Post-%s",mn_tag)),levels=c("Pre-model",sprintf("Post-%s",mn_tag))),
    I=c(moran_pre$estimate[1],moran_post$estimate[1]))
  mdf$col<-c("#D6604D","#2166AC");mdf$lbl<-sprintf("I=%+.4f",mdf$I)
  p3<-ggplot(mdf,aes(x=단계,y=I,fill=col))+geom_col(width=0.45,alpha=0.85)+
    geom_hline(yintercept=0,linetype="dashed")+geom_text(aes(label=lbl),vjust=-0.5,size=5)+
    scale_fill_identity()+labs(title="Moran's I — Spatial Autocorrelation",x="",y="Moran's I")+
    theme_minimal(base_family=FONT)
  ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_MoransI_%s.png",TS)),p3,width=7,height=6,dpi=150)
  ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_MoransI_%s.tiff",TS)),p3,width=7,height=6,dpi=300,compression="lzw")
  cat("    ✅ Moran's I (PNG+TIFF)\n")}

# PLOT 4a: 공간효과 지도
cat("  [Fig 4] 공간효과 + 고위험분류\n")
if(!is.null(bm)&&!is.null(bm$summary.random$idarea)){
  na<-nrow(shp_main);sp_eff<-bm$summary.random$idarea$mean[1:na]
  sp_lo2<-bm$summary.random$idarea$`0.025quant`[1:na];sp_hi2<-bm$summary.random$idarea$`0.975quant`[1:na]
  shp_plot<-shp_main;shp_plot$sp<-sp_eff;lim<-max(abs(sp_eff),na.rm=TRUE)
  p4a<-ggplot(shp_plot)+geom_sf(aes(fill=sp),color="white",linewidth=0.15)+
    scale_fill_gradient2(low="#2166AC",mid="white",high="#D6604D",midpoint=0,limits=c(-lim,lim))+
    labs(title=sprintf("Spatial Random Effect (%s)",mn_tag))+theme_void(base_family=FONT)
  ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_공간효과_%s.png",TS)),p4a,width=10,height=8,dpi=150)
  ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_공간효과_%s.tiff",TS)),p4a,width=10,height=8,dpi=300,compression="lzw")
  cat("    ✅ 공간효과 (PNG+TIFF)\n")

  # PLOT 4b: 고위험/저위험 분류 지도
  shp_plot$risk<-factor(dplyr::case_when(sp_lo2>0~"High",sp_hi2<0~"Low",TRUE~"Non-sig"),
    levels=c("High","Non-sig","Low"))
  p4b<-ggplot(shp_plot)+geom_sf(aes(fill=risk),color="white",linewidth=0.15)+
    scale_fill_manual(values=c("High"="#D6604D","Non-sig"="#F0F0F0","Low"="#2166AC"))+
    labs(title=sprintf("Risk Classification (%s)",mn_tag))+theme_void(base_family=FONT)
  ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_고위험분류_%s.png",TS)),p4b,width=10,height=8,dpi=150)
  ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_고위험분류_%s.tiff",TS)),p4b,width=10,height=8,dpi=300,compression="lzw")
  cat(sprintf("    ✅ 고위험분류 (High %d, Low %d) (PNG+TIFF)\n",sum(sp_lo2>0),sum(sp_hi2<0)))}

# PLOT 5: Forest Plot
cat("  [Fig 5] Forest plot\n")
if(!is.null(res_final$mv) && nrow(res_final$mv)>0){
  mc <- res_final$mv[res_final$mv$sig=="★" & res_final$mv$방향일치!="❌역방향",,drop=FALSE]
  if(nrow(mc)>0){
    ev <- ENG[mc$var_kr]; ev[is.na(ev)] <- mc$var_kr[is.na(ev)]; ev <- unname(ev)
    fp <- data.frame(var_en=ev, IRR=mc$IRR, lo=mc$lo, hi=mc$hi,
      dir=ifelse(mc$IRR>1,"Risk","Protective"), stringsAsFactors=FALSE)
    fp <- fp %>% arrange(IRR) %>% mutate(var_en=factor(var_en,levels=var_en))
    CI_CLIP <- 5; fp$hi_c <- pmin(fp$hi, CI_CLIP); fp$clipped <- fp$hi > CI_CLIP
    fp$irr_lab <- sprintf("%.2f (%.2f\u2013%.2f)", fp$IRR, fp$lo, fp$hi)
    p5 <- ggplot(fp, aes(x=IRR, y=var_en))+
      geom_vline(xintercept=1, linetype="21", color="gray50", linewidth=0.35)+
      geom_errorbarh(aes(xmin=lo, xmax=hi_c, color=dir), height=0, linewidth=0.7, show.legend=FALSE)+
      geom_point(aes(color=dir,fill=dir), shape=23, size=3, stroke=0.4)+
      geom_text(aes(label=irr_lab), hjust=-0.15, size=2.5, family=FONT, color="gray20")+
      scale_color_manual(name=NULL, values=c("Risk"=COL_RISK,"Protective"=COL_PROT),
        labels=c("Risk"="Risk (IRR > 1)","Protective"="Protective (IRR < 1)"))+
      scale_fill_manual(name=NULL, values=c("Risk"=COL_RISK,"Protective"=COL_PROT),
        labels=c("Risk"="Risk (IRR > 1)","Protective"="Protective (IRR < 1)"))+
      scale_x_log10()+
      labs(x="Incidence Rate Ratio (95% CrI)", y=NULL,
        title=sprintf("HAV — Significant Variables (%s, n=%d)",mn_tag,ifelse(is.null(res_final$N_final),0,res_final$N_final)))+
      theme_minimal(base_family=FONT, base_size=9)+
      theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        legend.position="bottom", plot.margin=margin(6,40,6,2,"mm"))
    ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_ForestPlot_%s.png",TS)),p5,width=190,height=max(120,30+nrow(fp)*14),units="mm",dpi=300,bg="white")
    ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_ForestPlot_%s.tiff",TS)),p5,width=190,height=max(120,30+nrow(fp)*14),units="mm",dpi=300,bg="white",compression="lzw")
    cat("    ✅ Forest plot (PNG+TIFF)\n")
  }
}

# PLOT 6: Spatiotemporal Interaction (a) 공간패턴 + (b) 시간추세
cat("  [Fig 6] Spatiotemporal interaction\n")
tryCatch({
  fit_final <- res_final$fit; ic_final <- res_final$ic
  if(!is.null(fit_final) && "idarea_time" %in% names(fit_final$summary.random)){
    psi <- fit_final$summary.random$idarea_time
    psi$idarea <- ic_final$idarea; psi$idtime <- ic_final$idtime
    psi$year <- ic_final$year; psi$region <- ic_final$region

    psi_spatial <- psi %>% group_by(region) %>%
      summarise(mean_psi=mean(mean, na.rm=TRUE), .groups="drop")
    shp_psi <- shp_main %>% left_join(psi_spatial, by="region")

    n_brk <- 5
    brk_psi <- quantile(shp_psi$mean_psi, probs=seq(0,1,length.out=n_brk+1), na.rm=TRUE)
    brk_psi <- unique(round(brk_psi, 3))
    if(length(brk_psi) < 3) brk_psi <- pretty(range(shp_psi$mean_psi, na.rm=TRUE), n=5)
    lbl_psi <- c()
    for(j in 1:(length(brk_psi)-1)) lbl_psi <- c(lbl_psi, sprintf("%.3f to %.3f", brk_psi[j], brk_psi[j+1]))
    shp_psi$psi_grp <- cut(shp_psi$mean_psi, breaks=brk_psi, labels=lbl_psi, include.lowest=TRUE)
    pal_psi <- colorRampPalette(c("#F7FBFF","#DEEBF7","#C6DBEF","#6BAED6","#2171B5","#08306B"))(length(lbl_psi))

    p6a <- ggplot(shp_psi) + geom_sf(aes(fill=psi_grp), color="gray45", linewidth=0.08) +
      scale_fill_manual(values=pal_psi, name=expression(paste("Mean ", psi)),
        drop=FALSE, na.value="gray88") +
      theme_void(base_family=FONT) + theme(legend.position=c(0.15,0.32))

    psi_temporal <- psi %>% group_by(year) %>%
      summarise(mean_psi=mean(mean, na.rm=TRUE),
        lo_psi=mean(`0.025quant`, na.rm=TRUE),
        hi_psi=mean(`0.975quant`, na.rm=TRUE), .groups="drop")

    p6b <- ggplot(psi_temporal, aes(x=year, y=mean_psi)) +
      geom_ribbon(aes(ymin=lo_psi, ymax=hi_psi), fill="gray80", alpha=0.6) +
      geom_line(linewidth=0.9, color="black") +
      geom_point(size=2, color="black") +
      geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=0.4) +
      scale_x_continuous(breaks=YEAR_START:YEAR_END) +
      labs(x="Year", y=expression(paste("Space-time interaction ", psi, " (posterior mean)"))) +
      theme_minimal(base_family=FONT, base_size=10) +
      theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=45, hjust=1))

    p6_combined <- tryCatch({
      p6_out <- p6a + p6b + plot_layout(widths=c(1.2, 1)) +
        plot_annotation(tag_levels="a", tag_prefix="(", tag_suffix=")")
      p6_out
    }, error=function(e){
      ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_Fig6a_STI_Spatial_%s.png",TS)),p6a,width=170,height=195,units="mm",dpi=300,bg="white")
      ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_Fig6b_STI_Temporal_%s.png",TS)),p6b,width=150,height=120,units="mm",dpi=300,bg="white")
      cat("    (patchwork 미설치 → 개별 저장)\n"); NULL
    })

    if(!is.null(p6_combined)){
      ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_STInteraction_%s.png",TS)),p6_combined,width=320,height=180,units="mm",dpi=300,bg="white")
      ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_STInteraction_%s.tiff",TS)),p6_combined,width=320,height=180,units="mm",dpi=300,bg="white",compression="lzw")
    }
    ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_Fig6a_STI_Spatial_%s.png",TS)),p6a,width=170,height=195,units="mm",dpi=300,bg="white")
    ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_Fig6a_STI_Spatial_%s.tiff",TS)),p6a,width=170,height=195,units="mm",dpi=300,bg="white",compression="lzw")
    ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_Fig6b_STI_Temporal_%s.png",TS)),p6b,width=150,height=120,units="mm",dpi=300,bg="white")
    ggsave(file.path(DIR_OUT,sprintf("HAV_v7.13_NODOG_Fig6b_STI_Temporal_%s.tiff",TS)),p6b,width=150,height=120,units="mm",dpi=300,bg="white",compression="lzw")
    cat("    ✅ STI (a) spatial + (b) temporal (PNG+TIFF)\n")
  } else { cat("    ⚠️ idarea_time 효과 없음 → STI 생략\n") }
}, error=function(e) cat(sprintf("    ❌ STI: %s\n", e$message)))

cat("\n  ═══ Figure 완료 ═══\n")
}, error=function(e) cat(sprintf("\n  ❌ Fig: %s\n", e$message)))


# ══════════════════════════════════════════
# PART 9. 잔차/프레딕션 지도
# ══════════════════════════════════════════
cat("\n## PART 9. 잔차/프레딕션 지도\n\n")
tryCatch({
dtag<-"HAV"; mn_tag<-sprintf("M%d",bi)
cat(sprintf("  질병: %s | 모델: %s\n",DISEASE_NAME,mn_tag))

obs_rate<-cor_merged%>%filter(region%in%shp_main$region,population>0)%>%
  group_by(region)%>%summarise(obs_cases=sum(cases,na.rm=TRUE),obs_pop=sum(population,na.rm=TRUE),
    obs_rate=obs_cases/obs_pop*100000,.groups="drop")

if(!is.null(bm)&&!is.null(bm$summary.fitted.values)){
  n_fit<-min(nrow(ic),nrow(bm$summary.fitted.values))
  zinb_correction <- 1  # NB → ZI 보정 없음 (곱하기 1)
  cat(sprintf("  ★ NB 예측: ZI 보정 없음 (직접 fitted values 사용)\n"))
  pred_df<-data.frame(region=ic$region[1:n_fit],population=ic$population[1:n_fit],
    fitted=bm$summary.fitted.values$mean[1:n_fit] * zinb_correction)
  pred_rate<-pred_df%>%group_by(region)%>%summarise(pred_cases=sum(fitted,na.rm=TRUE),
    pred_pop=sum(population,na.rm=TRUE),pred_rate=pred_cases/pred_pop*100000,.groups="drop")
  # ★ 예측값 캡핑 (v7.2: 99.5th percentile)
  cap_val <- quantile(pred_rate$pred_rate, 0.995, na.rm=TRUE)
  cap_val <- max(cap_val, max(obs_rate$obs_rate, na.rm=TRUE) * 2, 0.01)
  n_capped <- sum(pred_rate$pred_rate > cap_val, na.rm=TRUE)
  if(n_capped > 0) {
    cat(sprintf("  ⚠️ 예측값 캡핑: %d개 시군구 (cap=%.1f/100k)\n", n_capped, cap_val))
    for(j in which(pred_rate$pred_rate > cap_val)) {
      cat(sprintf("    %s: pred=%.1e → %.1f\n", pred_rate$region[j], pred_rate$pred_rate[j], cap_val))
    }
    pred_rate$pred_rate <- pmin(pred_rate$pred_rate, cap_val)
    pred_rate$pred_cases <- pred_rate$pred_rate * pred_rate$pred_pop / 100000
  }
  map_data<-obs_rate%>%left_join(pred_rate,by="region")%>%
    mutate(residual_rate=obs_rate-pred_rate,residual_ratio=ifelse(pred_rate>0,obs_rate/pred_rate,NA))
  shp_plot2<-shp_main%>%left_join(map_data,by="region")%>%
    mutate(obs_rate=ifelse(is.na(obs_rate),0,obs_rate),pred_rate=ifelse(is.na(pred_rate),0,pred_rate),
      residual_rate=ifelse(is.na(residual_rate),0,residual_rate))
  cat(sprintf("  시군구: %d | Obs: %.1f~%.1f | Pred: %.1f~%.1f | Resid: %+.1f~%+.1f\n",
    nrow(shp_plot2),min(shp_plot2$obs_rate),max(shp_plot2$obs_rate),
    min(shp_plot2$pred_rate),max(shp_plot2$pred_rate),
    min(shp_plot2$residual_rate),max(shp_plot2$residual_rate)))
  rate_max2<-max(quantile(shp_plot2$obs_rate,0.98,na.rm=TRUE),quantile(shp_plot2$pred_rate,0.98,na.rm=TRUE),0.01)

  p_obs<-ggplot(shp_plot2)+geom_sf(aes(fill=obs_rate),color="white",linewidth=0.15)+
    scale_fill_gradientn(colors=c("#FFF7EC","#FEE8C8","#FDD49E","#FDBB84","#FC8D59","#E34A33","#B30000"),
      limits=c(0,rate_max2),name="Rate\n/100k",oob=scales::squish)+
    labs(title="(A) Observed incidence rate")+theme_void(base_family=FONT)+
    theme(plot.title=element_text(size=13,face="bold",hjust=0.5))
  p_pred<-ggplot(shp_plot2)+geom_sf(aes(fill=pred_rate),color="white",linewidth=0.15)+
    scale_fill_gradientn(colors=c("#FFF7EC","#FEE8C8","#FDD49E","#FDBB84","#FC8D59","#E34A33","#B30000"),
      limits=c(0,rate_max2),name="Rate\n/100k",oob=scales::squish)+
    labs(title=sprintf("(B) Predicted rate (%s)",mn_tag))+theme_void(base_family=FONT)+
    theme(plot.title=element_text(size=13,face="bold",hjust=0.5))
  res_lim<-max(abs(quantile(shp_plot2$residual_rate,c(0.02,0.98),na.rm=TRUE)))
  p_res<-ggplot(shp_plot2)+geom_sf(aes(fill=residual_rate),color="white",linewidth=0.15)+
    scale_fill_gradient2(low="#2166AC",mid="#F7F7F7",high="#B2182B",midpoint=0,
      limits=c(-res_lim,res_lim),name="Residual\n/100k",oob=scales::squish)+
    labs(title="(C) Residual (observed - predicted)")+theme_void(base_family=FONT)+
    theme(plot.title=element_text(size=13,face="bold",hjust=0.5))

  p_combined<-p_obs+p_pred+p_res+plot_layout(nrow=1,guides="collect")+
    plot_annotation(title=sprintf("%s — Observed vs Predicted vs Residual (%s, %d-%d)",DISEASE_NAME,mn_tag,YEAR_START,YEAR_END),
      theme=theme(plot.title=element_text(size=15,face="bold",hjust=0.5)))
  fn_3p<-file.path(DIR_OUT,sprintf("%s_잔차프레딕션_3종_%s.png",dtag,TS))
  ggsave(fn_3p,p_combined,width=28,height=9,dpi=150)
  ggsave(file.path(DIR_OUT,sprintf("%s_Observed_%s.png",dtag,TS)),p_obs,width=10,height=8,dpi=150)
  ggsave(file.path(DIR_OUT,sprintf("%s_Predicted_%s_%s.png",dtag,mn_tag,TS)),p_pred,width=10,height=8,dpi=150)
  ggsave(file.path(DIR_OUT,sprintf("%s_Residual_%s_%s.png",dtag,mn_tag,TS)),p_res,width=10,height=8,dpi=150)
  cat(sprintf("  ✅ 3종 세트\n  ✅ Observed\n  ✅ Predicted\n  ✅ Residual\n"))

  # 연도별 비교
  cat("\n  연도별 비교 지도...\n")
  yearly_obs<-cor_merged%>%filter(region%in%shp_main$region,population>0)%>%mutate(rate=cases/population*100000)
  yearly_pred<-pred_df%>%mutate(year=ic$year[1:n_fit])%>%mutate(rate=fitted/population*100000)
  yr_max<-max(quantile(yearly_obs$rate,0.98,na.rm=TRUE),quantile(yearly_pred$rate,0.98,na.rm=TRUE),0.01)
  obs_list<-list();pred_list<-list()
  for(yr in years_use){
    yd_obs<-yearly_obs%>%filter(year==yr)%>%dplyr::select(region,rate)
    sy_obs<-shp_main%>%left_join(yd_obs,by="region")%>%mutate(rate=ifelse(is.na(rate),0,rate))
    obs_list[[as.character(yr)]]<-ggplot(sy_obs)+geom_sf(aes(fill=rate),color="white",linewidth=0.08)+
      scale_fill_gradientn(colors=c("#FFF7EC","#FDD49E","#FC8D59","#B30000"),limits=c(0,yr_max),name="/100k",oob=scales::squish)+
      labs(title=sprintf("%d Obs",yr))+theme_void(base_family=FONT)+theme(plot.title=element_text(size=10))
    yd_pred<-yearly_pred%>%filter(year==yr)%>%group_by(region)%>%summarise(rate=mean(rate,na.rm=TRUE),.groups="drop")
    sy_pred<-shp_main%>%left_join(yd_pred,by="region")%>%mutate(rate=ifelse(is.na(rate),0,rate))
    pred_list[[as.character(yr)]]<-ggplot(sy_pred)+geom_sf(aes(fill=rate),color="white",linewidth=0.08)+
      scale_fill_gradientn(colors=c("#FFF7EC","#FDD49E","#FC8D59","#B30000"),limits=c(0,yr_max),name="/100k",oob=scales::squish)+
      labs(title=sprintf("%d Pred",yr))+theme_void(base_family=FONT)+theme(plot.title=element_text(size=10))
  }
  ny<-length(years_use)
  p_yr<-wrap_plots(c(obs_list,pred_list),nrow=2,ncol=ny)+plot_layout(guides="collect")+
    plot_annotation(title=sprintf("%s — Observed (top) vs Predicted (bottom), %s",DISEASE_NAME,mn_tag),
      theme=theme(plot.title=element_text(size=14,face="bold")))
  ggsave(file.path(DIR_OUT,sprintf("%s_연도별_ObsPred_%s.png",dtag,TS)),p_yr,width=22,height=14,dpi=120)
  cat("  ✅ 연도별 Obs vs Pred\n")

  # 잔차 요약
  cat("\n  잔차 요약:\n")
  cat(sprintf("  Mean: %+.3f | Median: %+.3f | SD: %.3f\n",
    mean(shp_plot2$residual_rate,na.rm=TRUE),median(shp_plot2$residual_rate,na.rm=TRUE),sd(shp_plot2$residual_rate,na.rm=TRUE)))
  top_over<-map_data%>%arrange(desc(residual_rate))%>%head(5)
  top_under<-map_data%>%arrange(residual_rate)%>%head(5)
  cat("  과소예측 Top5:\n")
  for(j in 1:nrow(top_over))cat(sprintf("    %s: obs=%.1f pred=%.1f resid=%+.1f\n",top_over$region[j],top_over$obs_rate[j],top_over$pred_rate[j],top_over$residual_rate[j]))
  cat("  과대예측 Top5:\n")
  for(j in 1:nrow(top_under))cat(sprintf("    %s: obs=%.1f pred=%.1f resid=%+.1f\n",top_under$region[j],top_under$obs_rate[j],top_under$pred_rate[j],top_under$residual_rate[j]))
}else{cat("  ❌ 모델 없음\n")}
},error=function(e)cat(sprintf("  ❌ 잔차지도: %s\n",e$message)))


# ══════════════════════════════════════════
# 완료
# ══════════════════════════════════════════
cat(sprintf("\n✅ 로그: %s\n",LOG))
sink()
cat(sprintf("\n═══ HAV v7.12 NB + 통합회귀 + 가설강화 + AUTO + Academic Output 완료 ═══\n"))
cat(sprintf("★ 모델: %s (Negative Binomial)\n", FAMILY))
cat(sprintf("★ 목표: 정방향≥%d + 역방향=0 + Moran p>0.05\n", TARGET_FWD))
cat(sprintf("★ 최종: 정방향=%d 역방향=%d | 변수 %d개\n", n_fwd_final, n_rev_final, nrow(TV_FINAL)))
cat(sprintf("★ 엑셀: %s\n", if(exists("fn_xlsx")) fn_xlsx else ""))
cat(sprintf("★ Word: %s\n", if(exists("fn_docx")) fn_docx else ""))

# ═══════════════════════════════════════════════════════════════════════════════
# ▣ [통합] HAV 추가분석 — 동반 스크립트 4종을 같은 세션에서 이어 실행
#   메인(위 PART 1~9)이 만든 res_final / g_main / ic / shp_main 을 그대로 이어받음.
#   산출물은 메인과 같은 폴더(DIR_OUT)로.
# ═══════════════════════════════════════════════════════════════════════════════
DIR_LOG <- if (exists("DIR_OUT")) DIR_OUT else getwd()
OUTPUT_DIR <- DIR_LOG
if (!exists("res_final")) stop("메인(HAV v7.13_NODOG)을 먼저 끝까지 실행해야 합니다 (res_final 필요)")
cat("\n\n", strrep("█",80), "\n  통합 추가분석 시작 (8prior · φ계열 · 8graph · GiStar/LISA)\n", strrep("█",80), "\n", sep="")


#==============================================================================
# ▶ 통합블록: BYM2 PC-prior 8설정 민감도 = Table 2-3.S3
#==============================================================================
# ════════════════════════════════════════════════════════════════════════════════
# HAV — BYM2 PC-prior 민감도 분석 (종심 코멘트 1-3, 유대성 교수님)
# ════════════════════════════════════════════════════════════════════════════════
# 목적 : 본문 HAV BYM2 모형의 PC prior 가 지나치게 보수적이어서 φ(phi)≈1 이
#        prior 탓일 수 있다는 지적에 대해, prior 를 여러 설정으로 바꿔가며
#        (a) φ 사후분포, (b) DIC/WAIC, (c) 모든 공변량 IRR 의 변화를 표로 보고.
#        → "결과가 prior 설정에 둔감(robust)함" 또는 "민감함"을 객관적으로 입증.
#
# 사용법 : HAV 메인 분석(v7.13_NODOG / v8 AUTO)을 먼저 돌려
#          res_final, g_main, ic, shp_main 객체가 메모리에 있는 상태에서 source().
#          (= HAV_aux_BYM2.R 의 A1 블록과 동일한 사전조건)
#
# 현재 본문(reference) prior :
#          prec : pc.prec  param=c(1, 0.01)   → P(σ > 1)   = 0.01  (강한 수축)
#          phi  : pc       param=c(0.5, 2/3)  → P(φ < 0.5) = 2/3   (iid 선호)
#
# 산출 : HAV_PriorSensitivity_BYM2_<ts>.xlsx
#          00_prior_grid       — 시나리오 정의
#          01_phi_by_prior     — 시나리오별 φ 사후(mean·sd·2.5%·50%·97.5%)
#          02_fit_by_prior     — 시나리오별 DIC·WAIC·유의변수 수 (+ default 대비 Δ)
#          03_IRR_by_prior     — 공변량 × 시나리오 IRR (wide)
#          04_robustness       — default 대비 부호반전·유의성반전·최대 |Δ%| 요약
# Author : S.K.  /  Created : 2026-06-08
# ════════════════════════════════════════════════════════════════════════════════

tryCatch({
suppressMessages({
  library(INLA); library(dplyr); library(tidyr); library(openxlsx)
})

# ─── STEP 0. 사전조건 점검 ───────────────────────────────────────────────────
stopifnot(all(sapply(c("res_final","g_main","ic"), exists)))
ic <- res_final$ic
fixed_names <- rownames(res_final$fit$summary.fixed)
cov_main    <- fixed_names[fixed_names != "(Intercept)"]
cov_str     <- paste(paste0("`", cov_main, "`"), collapse = " + ")
ts <- format(Sys.time(), "%y%m%d_%H%M")
cat(sprintf("[HAV prior-sensitivity] covariates=%d  rows=%d  ts=%s\n",
            length(cov_main), nrow(ic), ts))

# 본문 A1 BYM2 와 동일하게 시간·상호작용 항 prior 를 고정(PC.prec(1,0.01)).
# 민감도 분석에서는 관심 하이퍼파라미터(BYM2 의 prec·phi)만 변화시키고 나머지는 본문 명세를 그대로 둔다.
# → 이래야 S1(principal) 행이 본문 BYM2 의 정규 DIC 를 재현하고 보충자료와 정합.
pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

# ─── STEP 1. prior 그리드 정의 ───────────────────────────────────────────────
# prec = c(U, alpha)  →  pc.prec : P(sigma > U) = alpha   (U↑ 또는 alpha↑ = 더 약한 수축)
# phi  = c(U, alpha)  →  pc      : P(phi   < U) = alpha   (alpha↓ = 공간구조 더 허용)
# "default" = INLA 내장 디폴트 사용(hyper 미지정)
prior_grid <- list(
  S1_default      = list(label="본문 reference (prec 1/0.01, phi 0.5/0.67)", prec=c(1,0.01),  phi=c(0.5, 2/3)),
  S2_prec_a05     = list(label="정밀도 prior 완화 (alpha 0.05)",            prec=c(1,0.05),  phi=c(0.5, 2/3)),
  S3_prec_U3      = list(label="정밀도 prior 완화 (U=3, alpha 0.05)",       prec=c(3,0.05),  phi=c(0.5, 2/3)),
  S4_prec_U5      = list(label="정밀도 prior 더 완화 (U=5, alpha 0.05)",    prec=c(5,0.05),  phi=c(0.5, 2/3)),
  S5_phi_neutral  = list(label="phi 중립 (P(phi<0.5)=0.5)",                 prec=c(1,0.01),  phi=c(0.5, 0.5)),
  S6_phi_spatial  = list(label="phi 공간 허용 (P(phi<0.5)=1/3)",            prec=c(1,0.01),  phi=c(0.5, 1/3)),
  S7_both_weak    = list(label="정밀도+phi 모두 완화",                       prec=c(3,0.05),  phi=c(0.5, 0.5)),
  S8_inla_default = list(label="INLA 내장 디폴트(hyper 미지정)",            prec="default",  phi="default")
)

build_hyper <- function(s) {
  if (identical(s$prec, "default")) return(NULL)        # NULL → hyper 미지정(디폴트)
  list(prec = list(prior="pc.prec", param=s$prec),
       phi  = list(prior="pc",      param=s$phi))
}

fit_one <- function(s) {
  hy <- build_hyper(s)
  f_bym2 <- if (is.null(hy)) {
    "+ f(idarea, model='bym2', graph=g_main, scale.model=TRUE, constr=TRUE)"
  } else {
    "+ f(idarea, model='bym2', graph=g_main, scale.model=TRUE, constr=TRUE, hyper=hy)"
  }
  form <- as.formula(paste(
    "cases ~", cov_str,
    "+ offset(log(population + 1))",
    f_bym2,
    "+ f(idtime, model='rw1', hyper=pc_prec)",
    "+ f(idarea_time, model='iid', hyper=pc_prec)"
  ))
  tryCatch(
    inla(form, family="nbinomial", data=ic,
         control.compute   = list(dic=TRUE, waic=TRUE, cpo=TRUE),
         control.predictor = list(link=1, compute=TRUE),
         control.inla      = list(strategy="adaptive", int.strategy="auto"),
         verbose=FALSE),
    error=function(e){ cat(sprintf("  ❌ %s: %s\n", s$label, e$message)); NULL })
}

# ─── STEP 2. 시나리오별 적합 ─────────────────────────────────────────────────
fits <- list(); phi_tab <- list(); fit_tab <- list(); irr_tab <- list()
for (nm in names(prior_grid)) {
  s <- prior_grid[[nm]]
  cat(sprintf("\n[%s] %s\n", nm, s$label))
  t0 <- Sys.time()
  fit <- fit_one(s)
  if (is.null(fit) || is.null(fit$dic$dic) || is.na(fit$dic$dic)) { cat("  skipped\n"); next }
  el <- as.numeric(difftime(Sys.time(), t0, units="mins"))
  fits[[nm]] <- fit

  # φ 사후
  phi_row <- grep("phi", rownames(fit$summary.hyperpar), ignore.case=TRUE, value=TRUE)
  if (length(phi_row) > 0) {
    pr <- fit$summary.hyperpar[phi_row[1], , drop=FALSE]
    phi_tab[[nm]] <- data.frame(scenario=nm, label=s$label,
                                phi_mean=round(pr[["mean"]],4),  phi_sd=round(pr[["sd"]],4),
                                phi_lo=round(pr[["0.025quant"]],4), phi_med=round(pr[["0.5quant"]],4),
                                phi_hi=round(pr[["0.975quant"]],4), stringsAsFactors=FALSE)
  }
  # 적합도
  fit_tab[[nm]] <- data.frame(scenario=nm, label=s$label,
                              DIC=round(fit$dic$dic,2), WAIC=round(fit$waic$waic,2),
                              minutes=round(el,2), stringsAsFactors=FALSE)
  # IRR
  fe <- fit$summary.fixed; fe <- fe[rownames(fe)!="(Intercept)", , drop=FALSE]
  irr_tab[[nm]] <- data.frame(
    variable = rownames(fe),
    IRR      = round(exp(fe[["mean"]]),4),
    lo       = round(exp(fe[["0.025quant"]]),4),
    hi       = round(exp(fe[["0.975quant"]]),4),
    sig      = ifelse(fe[["0.025quant"]]>0 | fe[["0.975quant"]]<0, "★",""),
    scenario = nm, stringsAsFactors=FALSE)
  cat(sprintf("  ✅ DIC=%.2f WAIC=%.2f  phi=%s  (%.1f min)\n",
              fit$dic$dic, fit$waic$waic,
              if (length(phi_row)>0) sprintf("%.3f",fit$summary.hyperpar[phi_row[1],"mean"]) else "NA", el))
}

phi_df <- bind_rows(phi_tab); fit_df <- bind_rows(fit_tab); irr_long <- bind_rows(irr_tab)
if (nrow(fit_df) > 0) { fit_df$dDIC <- fit_df$DIC - fit_df$DIC[1]; fit_df$dWAIC <- fit_df$WAIC - fit_df$WAIC[1] }

# IRR wide (변수 × 시나리오)
# 주의: select/mutate 는 INLA/MASS 와 충돌하므로 dplyr:: 명시 (NB→EHEC 디버그 lesson)
irr_wide <- irr_long %>%
  dplyr::mutate(cell = sprintf("%.3f (%.3f-%.3f)%s", IRR, lo, hi, ifelse(sig=="★"," ★",""))) %>%
  dplyr::select(variable, scenario, cell) %>%
  tidyr::pivot_wider(names_from=scenario, values_from=cell)

# ─── STEP 3. 강건성 요약 (default 대비) ─────────────────────────────────────
robust <- NULL
if (length(fits) > 1 && "S1_default" %in% names(fits)) {
  base <- irr_tab[["S1_default"]]
  rows <- list()
  for (nm in setdiff(names(irr_tab), "S1_default")) {
    cur <- irr_tab[[nm]]; m <- merge(base, cur, by="variable", suffixes=c("_base","_cur"))
    sign_flip <- sum((m$IRR_base>1) != (m$IRR_cur>1))
    sig_flip  <- sum(m$sig_base != m$sig_cur)
    maxd      <- max(abs(100*(m$IRR_cur-m$IRR_base)/m$IRR_base), na.rm=TRUE)
    rows[[nm]] <- data.frame(scenario=nm, label=prior_grid[[nm]]$label,
                             sign_flips=sign_flip, sig_flips=sig_flip,
                             max_abs_pct_change=round(maxd,2), stringsAsFactors=FALSE)
  }
  robust <- bind_rows(rows)
}

# ─── STEP 4. 저장 ────────────────────────────────────────────────────────────
grid_df <- bind_rows(lapply(names(prior_grid), function(nm){
  s<-prior_grid[[nm]]
  data.frame(scenario=nm, label=s$label,
             prec=ifelse(identical(s$prec,"default"),"default",paste(s$prec,collapse=", ")),
             phi =ifelse(identical(s$phi ,"default"),"default",paste(s$phi ,collapse=", ")),
             stringsAsFactors=FALSE)}))

wb <- createWorkbook()
addWorksheet(wb,"00_prior_grid");   writeData(wb,"00_prior_grid",grid_df)
addWorksheet(wb,"01_phi_by_prior"); writeData(wb,"01_phi_by_prior",phi_df)
addWorksheet(wb,"02_fit_by_prior"); writeData(wb,"02_fit_by_prior",fit_df)
addWorksheet(wb,"03_IRR_by_prior"); writeData(wb,"03_IRR_by_prior",irr_wide)
if (!is.null(robust)) { addWorksheet(wb,"04_robustness"); writeData(wb,"04_robustness",robust) }
outxlsx <- file.path(get0("OUTPUT_DIR", ifnotfound=getwd()),
                     sprintf("HAV_PriorSensitivity_BYM2_%s.xlsx", ts))
saveWorkbook(wb, outxlsx, overwrite=TRUE)
saveRDS(fits, file.path(get0("OUTPUT_DIR", ifnotfound=getwd()),
                        sprintf("HAV_PriorSensitivity_BYM2_fits_%s.rds", ts)))

cat("\n", strrep("═",78), "\n", sep="")
cat(sprintf("  저장: %s\n", outxlsx))
if (nrow(phi_df)>0) { cat("  φ 범위: ",
   sprintf("%.3f ~ %.3f (시나리오 간)\n", min(phi_df$phi_mean), max(phi_df$phi_mean))) }
if (!is.null(robust)) cat(sprintf("  최대 부호반전=%d · 최대 유의성반전=%d (default 대비)\n",
                                  max(robust$sign_flips), max(robust$sig_flips)))
cat(strrep("═",78), "\n")
}, error=function(e) cat(sprintf("\n\u26a0 [\ud1b5\ud569\ube14\ub85d \uac74\ub108\ub700: prior 8설정 S3] %s\n", conditionMessage(e))))
# ── 본문 해석 가이드 ──────────────────────────────────────────────────────────
# · φ 사후가 시나리오 전반에서 1 근처로 유지되면 → "φ≈1 은 prior 가 아니라 데이터가
#   강한 공간구조를 지지하기 때문"이라는 근거 (유대성 코멘트에 대한 직접 답).
# · DIC/WAIC 와 주요 IRR(상수도·지하수·주거지·재정자립도 등)이 부호·유의성 모두 유지되면
#   → "추론은 prior 설정에 둔감하다(robust)" 고 본문 Strengths/Sensitivity 에 1~2문장.
# · 만약 약한 prior 에서 φ 가 눈에 띄게 내려가거나 IRR 이 흔들리면 → 본문 수정 필요.
# ════════════════════════════════════════════════════════════════════════════════


#==============================================================================
# ▶ 통합블록: φ prior 계열(PC/Beta(1,1)/Beta(0.5,0.5)) = Table 2-3.S7 (E-1)
#==============================================================================
# ════════════════════════════════════════════════════════════════════════════════
# HAV — BYM2 혼합모수 φ 사전분포 "계열" 민감도 (유대성 교수님 이메일 E-1, 종심 전 필수)
# ════════════════════════════════════════════════════════════════════════════════
# 요청 원문(요약):
#   "Φ에 ①PC prior(현재) ②Beta(1,1) uniform ③Beta(0.5,0.5) Jeffreys 세 가지를 비교
#    적합하여 supplementary 1표로 보고. Φ posterior 가 0.85~1.0 에서 안정적이면 본문
#    표현 유지, prior 에 따라 0.4~1.0 으로 흔들리면 'near-pure structured spatial
#    variation' 표현을 톤 다운."
#
# ※ 기존 8-설정 민감도(HAV_PriorSensitivity_BYM2_260608.R)와 별개입니다.
#    그쪽은 precision·phi 의 PC prior 하이퍼파라미터를 흔든 것이고,
#    이 스크립트는 φ 의 prior '계열' 자체(PC vs Beta 두 종)를 바꿉니다.
#    precision prior 는 본문값 PC.prec(1,0.01) 로 3개 시나리오 모두 고정합니다.
#
# 사용법 : HAV 메인 분석을 먼저 돌려 res_final, g_main, ic 객체가 메모리에 있는 상태에서
#          source("HAV_PhiPriorFamily_BYM2_260609.R")  (= 기존 prior 민감도 스크립트와 동일 사전조건)
#
# 산출 : HAV_PhiPriorFamily_BYM2_<ts>.xlsx
#          01_phi_by_priorfamily  — prior 계열별 φ 사후(mean·sd·2.5%·50%·97.5%) + DIC·WAIC
#          02_IRR_by_priorfamily  — 공변량 × prior계열 IRR (부호·유의성 보존 확인용)
#          03_verdict             — φ 범위로 본문 표현 유지/톤다운 자동 판정
# Author : S.K.  /  Created : 2026-06-09
# ════════════════════════════════════════════════════════════════════════════════

tryCatch({
suppressMessages({ library(INLA); library(dplyr); library(tidyr); library(openxlsx) })

# ─── STEP 0. 사전조건 ────────────────────────────────────────────────────────────
stopifnot(all(sapply(c("res_final","g_main","ic"), exists)))
ic          <- res_final$ic
fixed_names <- rownames(res_final$fit$summary.fixed)
cov_main    <- fixed_names[fixed_names != "(Intercept)"]
cov_str     <- paste(paste0("`", cov_main, "`"), collapse = " + ")
ts <- format(Sys.time(), "%y%m%d_%H%M")
cat(sprintf("[HAV φ prior-family] covariates=%d  rows=%d  ts=%s\n",
            length(cov_main), nrow(ic), ts))

# 본문 명세 고정 항(시간·상호작용·정밀도 prior). φ prior 만 시나리오별로 교체.
pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

# ─── STEP 1. φ prior '계열' 3종 정의 ─────────────────────────────────────────────
# INLA bym2 의 φ 는 내부적으로 theta = logit(φ) 로 표현됨.  φ = 1/(1+exp(-theta)).
# Beta(a,b) prior on φ 를 theta 스케일 log-밀도로 옮기면(야코비안 dφ/dθ=φ(1-φ) 포함):
#   log π(θ) = a·log φ + b·log(1-φ) − lbeta(a,b)
#   · Beta(1,1)  : log φ + log(1-φ)                      (uniform)
#   · Beta(.5,.5): 0.5·log φ + 0.5·log(1-φ) − log(pi)    (Jeffreys)
phi_pc    <- list(prior = "pc", param = c(0.5, 2/3))   # ① 현재 본문(PC prior)
phi_unif  <- list(prior = "expression:
                     phi = 1/(1+exp(-theta));
                     log_dens = log(phi) + log(1-phi);
                     return(log_dens);")               # ② Beta(1,1) uniform
phi_jeff  <- list(prior = "expression:
                     phi = 1/(1+exp(-theta));
                     log_dens = 0.5*log(phi) + 0.5*log(1-phi) - log(pi);
                     return(log_dens);")               # ③ Beta(0.5,0.5) Jeffreys

prior_family <- list(
  P1_PC_current  = list(label = "① PC prior (current, P(φ<0.5)=2/3)", phi = phi_pc),
  P2_Beta11_unif = list(label = "② Beta(1,1) uniform on φ",            phi = phi_unif),
  P3_Beta05_jeff = list(label = "③ Beta(0.5,0.5) Jeffreys on φ",       phi = phi_jeff)
)

# ─── STEP 2. 적합 함수 (precision = PC.prec(1,0.01) 고정, φ prior 만 교체) ────────
fit_one <- function(s) {
  hy <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), phi = s$phi)
  form <- as.formula(paste(
    "cases ~", cov_str,
    "+ offset(log(population + 1))",
    "+ f(idarea, model='bym2', graph=g_main, scale.model=TRUE, constr=TRUE, hyper=hy)",
    "+ f(idtime, model='rw1', hyper=pc_prec)",
    "+ f(idarea_time, model='iid', hyper=pc_prec)"
  ))
  tryCatch(
    inla(form, family = "nbinomial", data = ic,
         control.compute   = list(dic = TRUE, waic = TRUE),
         control.predictor = list(link = 1, compute = TRUE),
         control.inla      = list(strategy = "adaptive", int.strategy = "auto"),
         verbose = FALSE),
    error = function(e) { cat(sprintf("  ❌ %s: %s\n", s$label, e$message)); NULL })
}

# ─── STEP 3. 시나리오별 적합 ─────────────────────────────────────────────────────
phi_tab <- list(); irr_tab <- list()
for (nm in names(prior_family)) {
  s <- prior_family[[nm]]
  cat(sprintf("\n[%s] %s\n", nm, s$label))
  fit <- fit_one(s)
  if (is.null(fit) || is.null(fit$dic$dic) || is.na(fit$dic$dic)) { cat("  skipped\n"); next }

  phi_row <- grep("phi", rownames(fit$summary.hyperpar), ignore.case = TRUE, value = TRUE)
  pr <- fit$summary.hyperpar[phi_row[1], , drop = FALSE]
  phi_tab[[nm]] <- data.frame(
    prior_family = nm, label = s$label,
    phi_mean = round(pr[["mean"]], 4),       phi_sd  = round(pr[["sd"]], 4),
    phi_lo   = round(pr[["0.025quant"]], 4), phi_med = round(pr[["0.5quant"]], 4),
    phi_hi   = round(pr[["0.975quant"]], 4),
    DIC = round(fit$dic$dic, 2), WAIC = round(fit$waic$waic, 2),
    stringsAsFactors = FALSE)

  fe <- fit$summary.fixed; fe <- fe[rownames(fe) != "(Intercept)", , drop = FALSE]
  irr_tab[[nm]] <- data.frame(
    variable = rownames(fe),
    cell = sprintf("%.3f (%.3f-%.3f)%s",
                   exp(fe[["mean"]]), exp(fe[["0.025quant"]]), exp(fe[["0.975quant"]]),
                   ifelse(fe[["0.025quant"]] > 0 | fe[["0.975quant"]] < 0, " ★", "")),
    prior_family = nm, stringsAsFactors = FALSE)
  cat(sprintf("  ✅ DIC=%.2f WAIC=%.2f  φ=%.3f (%.3f–%.3f)\n",
              fit$dic$dic, fit$waic$waic, pr[["mean"]], pr[["0.025quant"]], pr[["0.975quant"]]))
}
phi_df  <- bind_rows(phi_tab)
irr_wide <- bind_rows(irr_tab) %>% tidyr::pivot_wider(names_from = prior_family, values_from = cell)

# ─── STEP 4. 자동 판정 (본문 표현 유지 vs 톤다운) ────────────────────────────────
verdict <- NULL
if (nrow(phi_df) > 0) {
  rng_lo <- min(phi_df$phi_mean); rng_hi <- max(phi_df$phi_mean)
  stable <- (rng_lo >= 0.85)      # 세 prior 모두 φ̄ ≥ 0.85 이면 안정
  verdict <- data.frame(
    phi_mean_min = rng_lo, phi_mean_max = rng_hi,
    decision = if (stable)
      "STABLE → 본문 'near-pure structured spatial variation' 표현 유지"
    else
      "SENSITIVE → 본문 표현을 'predominantly structured'로 톤다운",
    stringsAsFactors = FALSE)
}

# ─── STEP 5. 저장 ────────────────────────────────────────────────────────────────
wb <- createWorkbook()
addWorksheet(wb, "01_phi_by_priorfamily"); writeData(wb, "01_phi_by_priorfamily", phi_df)
addWorksheet(wb, "02_IRR_by_priorfamily"); writeData(wb, "02_IRR_by_priorfamily", irr_wide)
if (!is.null(verdict)) { addWorksheet(wb, "03_verdict"); writeData(wb, "03_verdict", verdict) }
outxlsx <- file.path(get0("OUTPUT_DIR", ifnotfound = getwd()),
                     sprintf("HAV_PhiPriorFamily_BYM2_%s.xlsx", ts))
saveWorkbook(wb, outxlsx, overwrite = TRUE)

cat("\n", strrep("═", 78), "\n", sep = "")
cat(sprintf("  저장: %s\n", outxlsx))
if (nrow(phi_df) > 0) {
  cat(sprintf("  φ̄ 범위(prior 3계열): %.3f ~ %.3f\n", min(phi_df$phi_mean), max(phi_df$phi_mean)))
  cat("  →", verdict$decision, "\n")
}
cat(strrep("═", 78), "\n")
}, error=function(e) cat(sprintf("\n\u26a0 [\ud1b5\ud569\ube14\ub85d \uac74\ub108\ub700: φ prior 계열 S7] %s\n", conditionMessage(e))))
# ── 보고 가이드 ──────────────────────────────────────────────────────────────────
# · 01_phi_by_priorfamily 시트가 그대로 supplementary 표 1개(유대성 요청)입니다.
#   행=3 prior 계열, 열=φ(mean·sd·2.5%·50%·97.5%)·DIC·WAIC.
# · 03_verdict 가 STABLE 이면 본문 표현 유지, SENSITIVE 면 'predominantly structured'로
#   톤다운 → 그 결과만 알려주시면 supplementary 표 docx 와 본문 문구를 제가 마무리합니다.
# · (참고) 주모형은 BYM 합성곱이며 φ 는 BYM2 재매개화 진단값임을 본문/응답서에 이미 명시.
# ════════════════════════════════════════════════════════════════════════════════


#==============================================================================
# ▶ 통합블록: 8개 이웃그래프(Queen·Rook·KNN) 강건성 = Table 2-3.S2
#==============================================================================
# ════════════════════════════════════════════════════════════════════
# HAV 공간 분석 — Spatial Weight (Neighbor) 민감도 분석 [v3 — 장준수 정석]
# ════════════════════════════════════════════════════════════════════
# 목적: 학위논문 1차 심사 대비 — 천병철 교수 지시사항 (2026-04-30 미팅)
#
# 비교 대상 (장준수 박사 2025 Supp Table S2.1 / S4.1 형식 일치):
#   W1: Queen contiguity (현재 main 분석, 채택)
#   W2: Rook contiguity (모서리 공유만)
#   W3: KNN k=2
#   W4: KNN k=3
#   W5: KNN k=4
#   W6: KNN k=5
#   W7: KNN k=6
#   W8: KNN k=7
#
# 출력: M6 (BYM + RW1 + IID, NB) 결과 비교 테이블 + Moran's I 8종 비교
# 전제: 메인 분석 (HAV_v7.12_NB_*.R) 이 이미 실행되어 res_final 객체가 메모리에 존재
#
# v3 변경 (2026-05-03):
#   - W matrix 4종 → 8종 (장준수 박사 K=2~7 + Queen + Rook)
#   - Moran's I robust 계산 유지 (length matching + na.exclude + moran.mc fallback)
#
# 작성: 2026-05-02 (v1) / 2026-05-03 (v2 Moran fix) / 2026-05-03 (v3 K=2-7 확장)
# ════════════════════════════════════════════════════════════════════

tryCatch({
library(spdep); library(INLA); library(sf); library(dplyr)

cat("\n══════════════════════════════════════════\n")
cat("  HAV — Spatial Weight 민감도 분석 시작 (v3)\n")
cat("  비교: Queen + Rook + KNN k=2~7 (8종)\n")
cat("══════════════════════════════════════════\n\n")

# ─────────────────────────────────────────────────
# STEP 1. 8개 Neighbor 객체 생성
# ─────────────────────────────────────────────────

cat("[STEP 1] Neighbor 객체 생성 중...\n")

coords <- st_coordinates(st_centroid(st_geometry(shp_main)))
n_districts <- nrow(shp_main)

# Contiguity-based
nb_queen <- poly2nb(shp_main, snap=0.01, queen=TRUE)
nb_rook  <- poly2nb(shp_main, snap=0.01, queen=FALSE)

# KNN (k=2 to 7, 장준수 박사 형식)
nb_knn2 <- knn2nb(knearneigh(coords, k=2), sym=TRUE)
nb_knn3 <- knn2nb(knearneigh(coords, k=3), sym=TRUE)
nb_knn4 <- knn2nb(knearneigh(coords, k=4), sym=TRUE)
nb_knn5 <- knn2nb(knearneigh(coords, k=5), sym=TRUE)
nb_knn6 <- knn2nb(knearneigh(coords, k=6), sym=TRUE)
nb_knn7 <- knn2nb(knearneigh(coords, k=7), sym=TRUE)

nb_list <- list(
  Queen=nb_queen, Rook=nb_rook,
  KNN2=nb_knn2, KNN3=nb_knn3, KNN4=nb_knn4,
  KNN5=nb_knn5, KNN6=nb_knn6, KNN7=nb_knn7
)

nb_stats <- data.frame(
  W = names(nb_list),
  total_links = sapply(nb_list, function(x) sum(card(x))),
  mean_neighbors = sapply(nb_list, function(x) round(mean(card(x)),2)),
  median_neighbors = sapply(nb_list, function(x) median(card(x))),
  max_neighbors = sapply(nb_list, function(x) max(card(x))),
  isolated_count = sapply(nb_list, function(x) sum(card(x)==0))
)
print(nb_stats)
cat("\n")

# ─────────────────────────────────────────────────
# STEP 2. 각 Neighbor에 대해 INLA graph 생성
# ─────────────────────────────────────────────────

cat("[STEP 2] INLA graph 변환 중...\n")
graph_files <- list()
for(wname in names(nb_list)){
  fpath <- sprintf("/tmp/hav_sens_%s.graph", tolower(wname))
  nb2INLA(nb_list[[wname]], file=fpath)
  graph_files[[wname]] <- inla.read.graph(fpath)
  cat(sprintf("  ✅ %s\n", wname))
}
cat("\n")

# ─────────────────────────────────────────────────
# STEP 3. 메인 분석 결과의 covariate set 추출
# ─────────────────────────────────────────────────

cat("[STEP 3] 메인 분석 covariate set 확인 중...\n")

if(!exists("res_final") || is.null(res_final$ic)) {
  stop("res_final 객체가 없거나 비어 있습니다. 메인 스크립트를 먼저 실행하세요.\n")
}

ic_main <- res_final$ic
FMAP_main <- res_final$FMAP
cov_str_main <- paste(FMAP_main$safe, collapse=" + ")
cat(sprintf("  ✅ Main 분석 변수 %d개, N=%d, sig=%d, DIC=%.2f\n",
            nrow(FMAP_main), nrow(ic_main), res_final$sig_count, res_final$dic))
cat(sprintf("  [진단] ic_main idarea unique = %d, shp_main districts = %d\n\n",
            length(unique(ic_main$idarea)), n_districts))

# ─────────────────────────────────────────────────
# STEP 4. 각 W matrix로 M6 모델 재실행 (8회)
# ─────────────────────────────────────────────────

cat("[STEP 4] 8개 W matrix로 M6 (BYM+RW1+IID, NB) 재실행 중...\n")
cat("    예상 소요: 약 8~16분 (각 1~2분)\n\n")

pc_bym <- list(
  prec.unstruct = list(prior="pc.prec", param=c(0.5,0.01)),
  prec.spatial  = list(prior="pc.prec", param=c(0.5,0.01))
)
pc_prec <- list(prec=list(prior="pc.prec", param=c(0.5,0.01)))

formula_template <- function(graph_obj){
  as.formula(paste(
    "cases ~", cov_str_main,
    "+ offset(log(population+1))",
    "+ f(idarea, model='bym', graph=graph_obj, scale.model=TRUE, hyper=pc_bym)",
    "+ f(idtime, model='rw1', hyper=pc_prec)",
    "+ f(idarea_time, model='iid', hyper=pc_prec)"
  ))
}

# Robust Moran's I 계산 함수
compute_moran_post <- function(fit_w, ic_main, nb_obj_w, n_shp){
  resid_raw <- ic_main$cases - fit_w$summary.fitted.values$mean
  agg <- aggregate(resid_raw, by=list(idarea=ic_main$idarea),
                   FUN=mean, na.rm=TRUE)
  agg <- agg[order(agg$idarea), ]; names(agg)[2] <- "resid_mean"
  resid_vec <- rep(NA_real_, n_shp)
  resid_vec[agg$idarea] <- agg$resid_mean
  nb_listw <- nb2listw(nb_obj_w, style="W", zero.policy=TRUE)

  m_test <- tryCatch(
    moran.test(resid_vec, nb_listw, zero.policy=TRUE, na.action=na.exclude),
    error = function(e) NULL)
  if(!is.null(m_test)) return(list(
    I = as.numeric(m_test$estimate["Moran I statistic"]),
    p = as.numeric(m_test$p.value), method = "moran.test"))

  m_mc <- tryCatch(
    moran.mc(resid_vec, nb_listw, nsim=999, zero.policy=TRUE, na.action=na.exclude),
    error = function(e) NULL)
  if(!is.null(m_mc)) return(list(
    I = as.numeric(m_mc$statistic),
    p = as.numeric(m_mc$p.value), method = "moran.mc"))

  return(list(I = NA, p = NA, method = "failed"))
}

sens_results <- list()
t_start_all <- Sys.time()

for(wname in names(graph_files)){
  cat(sprintf("──────── [%s] ──────────\n", wname))
  t_start <- Sys.time()
  graph_obj <- graph_files[[wname]]

  fit_w <- tryCatch(
    inla(formula_template(graph_obj),
         family = "nbinomial", data = ic_main,
         control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE),
         control.predictor = list(link=1), verbose = FALSE),
    error = function(e){ cat(sprintf("  ❌ %s 실패: %s\n", wname, e$message)); NULL })

  if(is.null(fit_w) || is.na(fit_w$dic$dic)) {
    cat(sprintf("  ⚠ %s 모델 수렴 실패\n\n", wname)); next
  }

  fe <- fit_w$summary.fixed
  fe <- fe[rownames(fe) != "(Intercept)", , drop=FALSE]

  irr_df <- data.frame(
    var_eng = gsub("_z$", "", rownames(fe)),
    IRR = round(exp(fe$mean), 4),
    lo  = round(exp(fe$`0.025quant`), 4),
    hi  = round(exp(fe$`0.975quant`), 4),
    sig = ifelse(fe$`0.025quant` > 0 | fe$`0.975quant` < 0, "★", ""),
    stringsAsFactors = FALSE)

  n_sig <- sum(irr_df$sig == "★")
  moran_res <- compute_moran_post(fit_w, ic_main, nb_list[[wname]], n_districts)
  bym_eff <- fit_w$summary.random$idarea[1:n_distinct(ic_main$idarea), ]
  hi_risk <- sum(bym_eff$`0.025quant` > 0)
  lo_risk <- sum(bym_eff$`0.975quant` < 0)

  sens_results[[wname]] <- list(
    fit = fit_w, irr_df = irr_df, n_sig = n_sig,
    DIC = fit_w$dic$dic, WAIC = fit_w$waic$waic,
    moran_I_post = round(moran_res$I, 4),
    moran_p_post = round(moran_res$p, 4),
    moran_method = moran_res$method,
    high_risk = hi_risk, low_risk = lo_risk
  )

  elapsed <- as.numeric(difftime(Sys.time(), t_start, units="secs"))
  moran_disp <- ifelse(is.na(moran_res$I), "NA",
                       sprintf("%+.4f (p=%.3f)", moran_res$I, moran_res$p))
  cat(sprintf("  ✅ DIC=%.2f | WAIC=%.2f | Sig=%d/%d | Moran=%s | Hi=%d Lo=%d (%.0fs)\n\n",
              fit_w$dic$dic, fit_w$waic$waic, n_sig, nrow(irr_df),
              moran_disp, hi_risk, lo_risk, elapsed))
}

elapsed_all <- as.numeric(difftime(Sys.time(), t_start_all, units="mins"))
cat(sprintf("\n총 소요: %.1f분\n\n", elapsed_all))

# ─────────────────────────────────────────────────
# STEP 5. 비교 테이블
# ─────────────────────────────────────────────────

cat("══════════════════════════════════════════\n")
cat("  민감도 분석 결과 종합 (8종 W matrix)\n")
cat("══════════════════════════════════════════\n\n")

sens_summary <- data.frame(
  W = names(sens_results),
  total_links = sapply(names(sens_results), function(x) sum(card(nb_list[[x]]))),
  mean_neighbors = sapply(names(sens_results), function(x) round(mean(card(nb_list[[x]])), 2)),
  DIC = sapply(sens_results, function(x) round(x$DIC, 2)),
  WAIC = sapply(sens_results, function(x) round(x$WAIC, 2)),
  Sig_count = sapply(sens_results, function(x) x$n_sig),
  Moran_I_post = sapply(sens_results, function(x) x$moran_I_post),
  Moran_p_post = sapply(sens_results, function(x) x$moran_p_post),
  HighRisk = sapply(sens_results, function(x) x$high_risk),
  LowRisk = sapply(sens_results, function(x) x$low_risk)
)
sens_summary$dDIC_vs_Queen <- round(sens_summary$DIC - sens_summary$DIC[1], 2)
print(sens_summary)
cat("\n")

# Queen-sig IRR 비교
queen_sig <- sens_results$Queen$irr_df[sens_results$Queen$irr_df$sig == "★", ]
if(nrow(queen_sig) > 0){
  cat("[Queen-sig 변수 8-W 간 IRR 비교]\n\n")
  irr_compare <- queen_sig[, "var_eng", drop=FALSE]
  for(wname in names(sens_results)){
    df <- sens_results[[wname]]$irr_df
    df_sub <- df[match(queen_sig$var_eng, df$var_eng), ]
    irr_compare[[paste0(wname, "_IRR")]] <- df_sub$IRR
    irr_compare[[paste0(wname, "_sig")]] <- df_sub$sig
  }
  print(irr_compare)
  cat("\n")
}

# ─────────────────────────────────────────────────
# STEP 6. 엑셀 저장 (8 sheets per W + summary)
# ─────────────────────────────────────────────────

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "S1_W민감도요약_8종")
addWorksheet(wb, "S2_Queen_sig_IRR비교")

writeData(wb, "S1_W민감도요약_8종", sens_summary)
if(nrow(queen_sig) > 0) writeData(wb, "S2_Queen_sig_IRR비교", irr_compare)

# 각 W별 IRR 시트
for(wname in names(sens_results)){
  sheet_name <- paste0("IRR_", wname)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, sens_results[[wname]]$irr_df)
}

OUT_XLSX <- file.path(getwd(), sprintf("HAV_Sensitivity_W_v3_%s.xlsx",
                                       format(Sys.time(), "%y%m%d_%H%M")))
saveWorkbook(wb, OUT_XLSX, overwrite=TRUE)
cat(sprintf("✅ 엑셀 저장: %s\n\n", OUT_XLSX))

# ─────────────────────────────────────────────────
# STEP 7. MD 보고서 저장 (G.Downloads)
# ─────────────────────────────────────────────────

OUT_MD <- file.path(DIR_OUT, "HAV_Sensitivity_W_matrix.md")

md_lines <- c(
  "# HAV 공간 분석 — Spatial Weight Matrix 민감도 분석 (v3)",
  "",
  paste0("**작성일**: ", format(Sys.time(), "%Y-%m-%d %H:%M"), " (v3 — 장준수 정석 8종)"),
  "**근거**: 천병철 교수 4/30 미팅 + 장준수 박사 (2025) Supp Table S4.1 형식",
  "**비교 대상**: Queen + Rook + KNN k=2~7 (총 8종)",
  "",
  "---",
  "",
  "## Supplementary Table S2-1b. Neighbor 객체 통계 (8종 W matrix)",
  "",
  "| Spatial weight | Total links | Mean nbrs | Median | Max | Isolated |",
  "|:---|---:|---:|---:|---:|---:|"
)

for(i in 1:nrow(nb_stats)){
  s <- nb_stats[i, ]
  md_lines <- c(md_lines, sprintf("| %s | %d | %.2f | %d | %d | %d |",
    s$W, s$total_links, s$mean_neighbors, s$median_neighbors, s$max_neighbors, s$isolated_count))
}

md_lines <- c(md_lines, "", "---", "",
  "## Supplementary Table S2-1c. Sensitivity of HAV Bayesian model fit to spatial weight specification (8 W matrices)",
  "",
  "| Spatial weight | DIC | WAIC | ΔDIC vs Queen | Sig vars | Moran I post | p-value | High-risk | Low-risk |",
  "|:---|---:|---:|---:|:---:|:---:|:---:|:---:|:---:|")

for(i in 1:nrow(sens_summary)){
  s <- sens_summary[i, ]
  moran_str <- ifelse(is.na(s$Moran_I_post), "—", sprintf("%+.3f", s$Moran_I_post))
  pval_str  <- ifelse(is.na(s$Moran_p_post), "—", sprintf("%.3f", s$Moran_p_post))
  md_lines <- c(md_lines, sprintf("| %s | %.2f | %.2f | %+.2f | %d | %s | %s | %d | %d |",
    s$W, s$DIC, s$WAIC, s$dDIC_vs_Queen, s$Sig_count,
    moran_str, pval_str, s$HighRisk, s$LowRisk))
}

if(exists("irr_compare") && nrow(irr_compare) > 0){
  md_lines <- c(md_lines, "", "---", "",
    "## Supplementary Table S2-1d. IRR estimates of Queen-significant covariates across 8 spatial weight specifications",
    "",
    paste0("| Variable |",
           paste(sapply(names(sens_results), function(w) sprintf(" %s IRR | sig |", w)), collapse="")),
    paste0("|:---|", paste(rep(":---:|:---:|", length(sens_results)), collapse="")))

  for(i in 1:nrow(irr_compare)){
    r <- irr_compare[i, ]
    row_str <- sprintf("| %s |", r$var_eng)
    for(wname in names(sens_results)){
      row_str <- paste0(row_str,
                       sprintf(" %.4f | %s |",
                               r[[paste0(wname, "_IRR")]],
                               r[[paste0(wname, "_sig")]]))
    }
    md_lines <- c(md_lines, row_str)
  }
}

md_lines <- c(md_lines, "", "---", "",
  "## Source",
  "- **분석 코드**: HAV_Sensitivity_W_matrix.R (v3 — 장준수 정석 8종)",
  paste0("- **메인 분석**: ", sprintf("%d covariates, N=%d, DIC=%.2f", nrow(FMAP_main), nrow(ic_main), res_final$dic)),
  "- **장준수 박사 reference**: Supplementary Table S2.1 / S4.1 (KNN k=2~7 + Queen)",
  "",
  "*End of HAV Spatial Weight Sensitivity Analysis (v3)*")

writeLines(md_lines, OUT_MD)
cat(sprintf("✅ MD 보고서 저장: %s\n", OUT_MD))

cat("\n✅ 민감도 분석 완료 (v3, 8종 W matrix)\n\n")
}, error=function(e) cat(sprintf("\n\u26a0 [\ud1b5\ud569\ube14\ub85d \uac74\ub108\ub700: 8graph W민감도 S2] %s\n", conditionMessage(e))))


#==============================================================================
# ▶ 통합블록: Getis-Ord Gi* 핫/콜드 + LISA = Figure 2-3.S2
#==============================================================================
# ════════════════════════════════════════════════════════════════════
# HAV 공간 분석 — Getis-Ord Gi* + Anselin LISA cluster (v1)
# ════════════════════════════════════════════════════════════════════
# 목적: 학위논문 1차 심사 대비 — 박사선배 4명 + 외부 3편 표준 정렬
# 권고안: HAV_권고안_P4.md ★★★ Must-add (Gi*) + ★★ Should-add (LISA, S2-2hi)
# 박사선배 표준: 장준수 Fig 2.4-2.5/3.3, 김지현 Fig 2-3
# 외부 표준: Boyce 2021 (EID) Fig 6, Rotondo 2018 (Sci Total Environ) Fig 3-4
#
# 비교 대상:
#   (A) Getis-Ord Gi* hotspot (z > 1.65/1.96/2.58 + 음의 cold spot)
#   (B) Anselin's LISA cluster (HH/LL/HL/LH at p<0.05)
#       (a) raw incidence rate
#       (b) covariate-adjusted BYM residual (M6 fit 결과 사용)
#
# 사전 조건 (메모리 보유):
#   res_final = run_model() 결과 list ($fit, $ic, $FMAP, ...)
#   shp_main  = sf object (223 districts, islands excluded)
#   nb_w      = nb2listw(nb_obj, style="W")
#
# 본 스크립트는 메인 분석 (HAV_v7_AUTO 또는 v8_AUTO_260418) 후 메모리 상태에서 실행됩니다.
# Sensitivity v3 (W matrix 8종) 와 동일한 객체 재사용 패턴.
#
# 작성: 2026-05-04
# ════════════════════════════════════════════════════════════════════

tryCatch({
suppressMessages({
  library(spdep); library(sf); library(dplyr); library(ggplot2)
  library(openxlsx); library(patchwork)
})

cat("\n", strrep("═", 78), "\n", sep = "")
cat("  HAV — Getis-Ord Gi* + Anselin LISA cluster analysis\n")
cat(strrep("═", 78), "\n\n", sep = "")

# ─────────────────────────────────────────────────
# STEP 0. Pre-flight 점검 (Sensitivity v3 와 동일 패턴)
# ─────────────────────────────────────────────────
cat("[STEP 0] Pre-flight 점검...\n")

required_objs <- c("res_final", "shp_main")
missing_objs  <- required_objs[!sapply(required_objs, exists)]
if (length(missing_objs) > 0)
  stop(sprintf("필수 객체 누락: %s\n메인 분석 (HAV_v7/v8 AUTO) 먼저 실행 필요",
               paste(missing_objs, collapse = ", ")))

if (is.null(res_final$ic) || nrow(res_final$ic) == 0) stop("res_final$ic 비어있음")
if (is.null(res_final$fit)) stop("res_final$fit 비어있음")

ic <- res_final$ic
n_districts <- nrow(shp_main)

# nb_w 가 메모리에 없으면 재생성
if (!exists("nb_w")) {
  nb_obj <- poly2nb(shp_main, snap = 0.01, queen = TRUE)
  nb_w   <- nb2listw(nb_obj, style = "W", zero.policy = TRUE)
  cat(sprintf("  ✅ nb_w 재생성 (Queen contiguity)\n"))
} else {
  cat(sprintf("  ✅ nb_w 재사용 (메모리)\n"))
}

cat(sprintf("  ✅ res_final: ic %d rows, fit DIC=%.2f\n", nrow(ic), res_final$fit$dic$dic))
cat(sprintf("  ✅ shp_main: %d districts (islands excluded)\n\n", n_districts))

# 단일 timestamp + 출력 폴더
TS <- format(Sys.time(), "%y%m%d_%H%M")
DIR_OUT <- if (exists("DIR_OUT_USER")) DIR_OUT_USER else
  file.path(getwd(), "output")
if (!dir.exists(DIR_OUT)) {
  tryCatch(dir.create(DIR_OUT, recursive = TRUE), error = function(e) {})
  if (!dir.exists(DIR_OUT))
    DIR_OUT <- file.path(Sys.getenv("HOME"), "Desktop")
}
cat(sprintf("  📁 DIR_OUT: %s\n", DIR_OUT))
cat(sprintf("  📅 TS: %s\n\n", TS))

# ─────────────────────────────────────────────────
# STEP 1. 5-year aggregate incidence (per 100,000) 계산
# ─────────────────────────────────────────────────
cat("[STEP 1] 5-year aggregate HAV incidence 계산...\n")

# ic 에서 region 별 aggregate
agg_inc <- ic %>%
  group_by(idarea) %>%
  summarise(
    cases       = sum(cases, na.rm = TRUE),
    pop_5y_mean = mean(population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(incidence_per_100k = cases / (pop_5y_mean * 5) * 100000) %>%
  arrange(idarea)

# shp_main 순서에 맞춰 vector 만들기 (idarea = 1..n_districts)
inc_vec <- rep(NA_real_, n_districts)
inc_vec[agg_inc$idarea] <- agg_inc$incidence_per_100k
inc_vec[is.na(inc_vec)] <- 0   # missing → 0 처리

shp_main$cases_5y    <- 0; shp_main$cases_5y[agg_inc$idarea]    <- agg_inc$cases
shp_main$pop_5y_mean <- 0; shp_main$pop_5y_mean[agg_inc$idarea] <- agg_inc$pop_5y_mean
shp_main$inc_5y      <- inc_vec

cat(sprintf("  📊 Incidence range: %.2f ~ %.2f / 100k (mean %.2f, SD %.2f)\n",
            min(inc_vec), max(inc_vec), mean(inc_vec), sd(inc_vec)))
cat(sprintf("  📊 N districts: %d (zero cases: %d)\n\n",
            n_districts, sum(inc_vec == 0)))

# ─────────────────────────────────────────────────
# STEP 2. Getis-Ord Gi* hotspot statistic
# ─────────────────────────────────────────────────
cat("[STEP 2] Getis-Ord Gi* hotspot 계산...\n")

# Gi* 는 self 포함 binary weight 사용
nb_obj_gi  <- if (exists("nb_obj")) nb_obj else poly2nb(shp_main, snap=0.01, queen=TRUE)
nb_self    <- include.self(nb_obj_gi)
W_self_B   <- nb2listw(nb_self, style = "B", zero.policy = TRUE)
gi_star    <- localG(inc_vec, W_self_B, zero.policy = TRUE)

shp_main$Gi_z <- as.numeric(gi_star)
shp_main$Gi_class <- factor(
  dplyr::case_when(
    shp_main$Gi_z >  2.58 ~ "Hot 99%",
    shp_main$Gi_z >  1.96 ~ "Hot 95%",
    shp_main$Gi_z >  1.65 ~ "Hot 90%",
    shp_main$Gi_z < -2.58 ~ "Cold 99%",
    shp_main$Gi_z < -1.96 ~ "Cold 95%",
    shp_main$Gi_z < -1.65 ~ "Cold 90%",
    TRUE                  ~ "Non-significant"
  ),
  levels = c("Cold 99%","Cold 95%","Cold 90%","Non-significant",
             "Hot 90%","Hot 95%","Hot 99%")
)

gi_summary <- table(shp_main$Gi_class)
cat("  Gi* class distribution:\n")
print(gi_summary)
cat(sprintf("  → Hot total: %d | Cold total: %d | NS: %d\n\n",
            sum(shp_main$Gi_class %in% c("Hot 90%","Hot 95%","Hot 99%")),
            sum(shp_main$Gi_class %in% c("Cold 90%","Cold 95%","Cold 99%")),
            sum(shp_main$Gi_class == "Non-significant")))

# ─────────────────────────────────────────────────
# STEP 3. Gi* 지도 출력 (PNG + TIFF)
# ─────────────────────────────────────────────────
cat("[STEP 3] Gi* hotspot map 생성...\n")

gi_pal <- c(
  "Cold 99%"="#08306b","Cold 95%"="#2171b5","Cold 90%"="#9ecae1",
  "Non-significant"="#f0f0f0",
  "Hot 90%"="#fcae91","Hot 95%"="#de2d26","Hot 99%"="#67000d"
)

p_gi <- ggplot(shp_main) +
  geom_sf(aes(fill = Gi_class), colour = "grey70", linewidth = 0.08) +
  scale_fill_manual(values = gi_pal, drop = FALSE,
                    name = "Gi* z-score\nclassification") +
  labs(
    title    = "Getis-Ord Gi* hotspot/coldspot of HAV foodborne disease incidence",
    subtitle = "South Korea, 2020-2024 (5-year aggregate, per 100,000)",
    caption  = "z thresholds: |z|=1.65 (p=0.10), 1.96 (p=0.05), 2.58 (p=0.01)"
  ) +
  theme_void(base_size = 10) +
  theme(
    legend.position = "right",
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, colour = "grey30"),
    plot.caption  = element_text(size = 8, colour = "grey50")
  )

fn_gi_png  <- file.path(DIR_OUT, sprintf("HAV_S2-3t_GiStar_HotspotMap_%s.png", TS))
fn_gi_tiff <- file.path(DIR_OUT, sprintf("HAV_S2-3t_GiStar_HotspotMap_%s.tiff", TS))
ggsave(fn_gi_png,  p_gi, width = 7, height = 8, dpi = 300, bg = "white")
ggsave(fn_gi_tiff, p_gi, width = 7, height = 8, dpi = 300, bg = "white", compression = "lzw")
cat(sprintf("  ✅ %s\n", basename(fn_gi_png)))
cat(sprintf("  ✅ %s\n\n", basename(fn_gi_tiff)))

# ─────────────────────────────────────────────────
# STEP 4. LISA panel (a): raw incidence rate
# ─────────────────────────────────────────────────
cat("[STEP 4] LISA cluster (a) raw incidence rate 계산...\n")

lm_raw  <- localmoran(inc_vec, nb_w, zero.policy = TRUE)
lag_raw <- lag.listw(nb_w, inc_vec, zero.policy = TRUE)
mean_inc <- mean(inc_vec, na.rm = TRUE)
mean_lag <- mean(lag_raw, na.rm = TRUE)

shp_main$LM_p_raw <- lm_raw[, "Pr(z != E(Ii))"]
shp_main$LISA_class_raw <- factor(
  dplyr::case_when(
    shp_main$LM_p_raw > 0.05                                    ~ "Non-significant",
    inc_vec > mean_inc & lag_raw > mean_lag                     ~ "High-High",
    inc_vec < mean_inc & lag_raw < mean_lag                     ~ "Low-Low",
    inc_vec > mean_inc & lag_raw < mean_lag                     ~ "High-Low",
    inc_vec < mean_inc & lag_raw > mean_lag                     ~ "Low-High",
    TRUE                                                        ~ "Non-significant"
  ),
  levels = c("High-High","Low-Low","High-Low","Low-High","Non-significant")
)

lisa_raw_summary <- table(shp_main$LISA_class_raw)
cat("  LISA (a) raw class distribution:\n")
print(lisa_raw_summary)
cat("\n")

# ─────────────────────────────────────────────────
# STEP 5. LISA panel (b): covariate-adjusted BYM residual (M6 fit)
# ─────────────────────────────────────────────────
cat("[STEP 5] LISA cluster (b) covariate-adjusted BYM residual 계산...\n")

# M6 fit 에서 fitted value 추출 → Pearson residual → district-mean
n_fit <- min(nrow(ic), nrow(res_final$fit$summary.fitted.values))
fitted_mean <- res_final$fit$summary.fitted.values$mean[1:n_fit]
pearson_res <- (ic$cases[1:n_fit] - fitted_mean) / sqrt(fitted_mean + 1e-6)

# district-level mean residual
res_df <- data.frame(idarea = ic$idarea[1:n_fit], r = pearson_res) %>%
  group_by(idarea) %>%
  summarise(r_mean = mean(r, na.rm = TRUE), .groups = "drop") %>%
  arrange(idarea)

resid_vec <- rep(NA_real_, n_districts)
resid_vec[res_df$idarea] <- res_df$r_mean
resid_vec[is.na(resid_vec)] <- 0

shp_main$residual_5y <- resid_vec
cat(sprintf("  Pearson residual range: %+.3f ~ %+.3f (mean %+.3f)\n",
            min(resid_vec), max(resid_vec), mean(resid_vec)))

lm_adj  <- localmoran(resid_vec, nb_w, zero.policy = TRUE)
lag_adj <- lag.listw(nb_w, resid_vec, zero.policy = TRUE)

shp_main$LM_p_adj <- lm_adj[, "Pr(z != E(Ii))"]
shp_main$LISA_class_adj <- factor(
  dplyr::case_when(
    shp_main$LM_p_adj > 0.05                  ~ "Non-significant",
    resid_vec > 0 & lag_adj > 0               ~ "High-High",
    resid_vec < 0 & lag_adj < 0               ~ "Low-Low",
    resid_vec > 0 & lag_adj < 0               ~ "High-Low",
    resid_vec < 0 & lag_adj > 0               ~ "Low-High",
    TRUE                                      ~ "Non-significant"
  ),
  levels = c("High-High","Low-Low","High-Low","Low-High","Non-significant")
)

lisa_adj_summary <- table(shp_main$LISA_class_adj)
cat("  LISA (b) BYM-residual class distribution:\n")
print(lisa_adj_summary)
cat("\n")

# ─────────────────────────────────────────────────
# STEP 6. LISA 2-panel 지도
# ─────────────────────────────────────────────────
cat("[STEP 6] LISA 2-panel 지도 생성...\n")

lisa_pal <- c(
  "High-High"="#de2d26","Low-Low"="#3182bd",
  "High-Low"="#fee08b","Low-High"="#bcbddc",
  "Non-significant"="#f0f0f0"
)

p_lisa_a <- ggplot(shp_main) +
  geom_sf(aes(fill = LISA_class_raw), colour = "grey70", linewidth = 0.08) +
  scale_fill_manual(values = lisa_pal, drop = FALSE, name = "LISA cluster") +
  labs(title = "(a) Raw HAV incidence rate") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10))

p_lisa_b <- ggplot(shp_main) +
  geom_sf(aes(fill = LISA_class_adj), colour = "grey70", linewidth = 0.08) +
  scale_fill_manual(values = lisa_pal, drop = FALSE, name = "LISA cluster") +
  labs(title = "(b) Covariate-adjusted BYM residual") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10))

p_lisa <- p_lisa_a + p_lisa_b +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Anselin's local Moran's I cluster classification",
    subtitle = "South Korea HAV foodborne disease, 2020-2024",
    theme    = theme(plot.title = element_text(face = "bold", size = 12))
  )

fn_lisa_png  <- file.path(DIR_OUT, sprintf("HAV_S2-3hi_LISA_2panel_%s.png", TS))
fn_lisa_tiff <- file.path(DIR_OUT, sprintf("HAV_S2-3hi_LISA_2panel_%s.tiff", TS))
ggsave(fn_lisa_png,  p_lisa, width = 12, height = 7, dpi = 300, bg = "white")
ggsave(fn_lisa_tiff, p_lisa, width = 12, height = 7, dpi = 300, bg = "white", compression = "lzw")
cat(sprintf("  ✅ %s\n", basename(fn_lisa_png)))
cat(sprintf("  ✅ %s\n\n", basename(fn_lisa_tiff)))

# ─────────────────────────────────────────────────
# STEP 7. xlsx export — district-level Gi* + LISA
# ─────────────────────────────────────────────────
cat("[STEP 7] xlsx Source data export...\n")

region_col <- if ("region" %in% names(shp_main)) shp_main$region else
              if ("SIG_KOR_NM" %in% names(shp_main)) shp_main$SIG_KOR_NM else
              paste0("idarea_", seq_len(n_districts))

out_tbl <- data.frame(
  idarea          = seq_len(n_districts),
  region          = region_col,
  cases_5y        = shp_main$cases_5y,
  pop_5y_mean     = round(shp_main$pop_5y_mean, 0),
  incidence_per_100k = round(shp_main$inc_5y, 3),
  Gi_z            = round(shp_main$Gi_z, 3),
  Gi_class        = as.character(shp_main$Gi_class),
  residual_5y     = round(shp_main$residual_5y, 3),
  LISA_class_raw  = as.character(shp_main$LISA_class_raw),
  LM_p_raw        = round(shp_main$LM_p_raw, 4),
  LISA_class_adj  = as.character(shp_main$LISA_class_adj),
  LM_p_adj        = round(shp_main$LM_p_adj, 4),
  stringsAsFactors = FALSE
)

summary_df <- data.frame(
  Metric = c("Hot 99% (Gi*)","Hot 95% (Gi*)","Hot 90% (Gi*)",
             "Cold 99% (Gi*)","Cold 95% (Gi*)","Cold 90% (Gi*)",
             "Non-significant (Gi*)",
             "LISA HH (raw)","LISA LL (raw)","LISA HL (raw)","LISA LH (raw)",
             "LISA HH (adj)","LISA LL (adj)","LISA HL (adj)","LISA LH (adj)"),
  Count  = c(sum(shp_main$Gi_class == "Hot 99%"),
             sum(shp_main$Gi_class == "Hot 95%"),
             sum(shp_main$Gi_class == "Hot 90%"),
             sum(shp_main$Gi_class == "Cold 99%"),
             sum(shp_main$Gi_class == "Cold 95%"),
             sum(shp_main$Gi_class == "Cold 90%"),
             sum(shp_main$Gi_class == "Non-significant"),
             sum(shp_main$LISA_class_raw == "High-High"),
             sum(shp_main$LISA_class_raw == "Low-Low"),
             sum(shp_main$LISA_class_raw == "High-Low"),
             sum(shp_main$LISA_class_raw == "Low-High"),
             sum(shp_main$LISA_class_adj == "High-High"),
             sum(shp_main$LISA_class_adj == "Low-Low"),
             sum(shp_main$LISA_class_adj == "High-Low"),
             sum(shp_main$LISA_class_adj == "Low-High"))
)

wb <- createWorkbook()
addWorksheet(wb, "01_Gi_LISA_district"); writeData(wb, "01_Gi_LISA_district", out_tbl)
addWorksheet(wb, "02_Summary");           writeData(wb, "02_Summary",         summary_df)

fn_xlsx <- file.path(DIR_OUT, sprintf("HAV_GiStar_LISA_summary_%s.xlsx", TS))
saveWorkbook(wb, fn_xlsx, overwrite = TRUE)
cat(sprintf("  ✅ %s\n\n", basename(fn_xlsx)))

# ─────────────────────────────────────────────────
# STEP 8. MD 보고서 (장준수 / Sensitivity v3 형식)
# ─────────────────────────────────────────────────
cat("[STEP 8] MD 보고서 생성...\n")

DIR_LOG <- "output"
OUT_MD <- file.path(DIR_LOG, "HAV_GiStar_LISA.md")

md <- c(
  "# HAV 공간 분석 — Getis-Ord Gi* + Anselin LISA cluster analysis",
  "",
  paste0("**작성일**: ", format(Sys.time(), "%Y-%m-%d %H:%M"), " (HAV_v7_AUTO 결과 기반)"),
  "**근거**: 박사선배 4명 표준 + 외부 3편 (Boyce 2021, Rotondo 2018, Mulder 2020)",
  "**대응**: HAV_권고안_P4.md ★★★ Must-add (Gi*) + ★★ Should-add (LISA)",
  "",
  "---",
  "",
  "## Methods",
  sprintf("- 분석 단위: %d 시군구 (도서 지역 제외)", n_districts),
  "- Outcome: 5-year aggregate HAV foodborne disease incidence per 100,000 person-years",
  "- Weight matrix: Queen contiguity (binary scheme for Gi*, row-standardised for LISA)",
  "- Gi* statistic: localG() with include.self() — z-score classification at |z|=1.65/1.96/2.58",
  "- LISA: localmoran() — High-High / Low-Low / High-Low / Low-High at p<0.05",
  "- 2-panel LISA: (a) raw incidence rate vs (b) M6 BYM Pearson residual (covariate-adjusted)",
  "",
  "---",
  "",
  "## Supplementary Figure S2-2t. Getis-Ord Gi* hotspot summary",
  "",
  "| Class | N districts | % |",
  "|:---|---:|---:|",
  sprintf("| Hot 99%% (z>2.58) | %d | %.1f%% |", sum(shp_main$Gi_class=="Hot 99%"),  100*sum(shp_main$Gi_class=="Hot 99%")/n_districts),
  sprintf("| Hot 95%% (z>1.96) | %d | %.1f%% |", sum(shp_main$Gi_class=="Hot 95%"),  100*sum(shp_main$Gi_class=="Hot 95%")/n_districts),
  sprintf("| Hot 90%% (z>1.65) | %d | %.1f%% |", sum(shp_main$Gi_class=="Hot 90%"),  100*sum(shp_main$Gi_class=="Hot 90%")/n_districts),
  sprintf("| Cold 99%% (z<-2.58) | %d | %.1f%% |", sum(shp_main$Gi_class=="Cold 99%"), 100*sum(shp_main$Gi_class=="Cold 99%")/n_districts),
  sprintf("| Cold 95%% (z<-1.96) | %d | %.1f%% |", sum(shp_main$Gi_class=="Cold 95%"), 100*sum(shp_main$Gi_class=="Cold 95%")/n_districts),
  sprintf("| Cold 90%% (z<-1.65) | %d | %.1f%% |", sum(shp_main$Gi_class=="Cold 90%"), 100*sum(shp_main$Gi_class=="Cold 90%")/n_districts),
  sprintf("| Non-significant | %d | %.1f%% |",   sum(shp_main$Gi_class=="Non-significant"), 100*sum(shp_main$Gi_class=="Non-significant")/n_districts),
  "",
  "---",
  "",
  "## Supplementary Figure S2-2hi. LISA cluster (a) raw vs (b) adjusted",
  "",
  "| LISA class | (a) Raw | (b) BYM-residual |",
  "|:---|---:|---:|",
  sprintf("| High-High | %d | %d |", sum(shp_main$LISA_class_raw=="High-High"), sum(shp_main$LISA_class_adj=="High-High")),
  sprintf("| Low-Low   | %d | %d |", sum(shp_main$LISA_class_raw=="Low-Low"),   sum(shp_main$LISA_class_adj=="Low-Low")),
  sprintf("| High-Low  | %d | %d |", sum(shp_main$LISA_class_raw=="High-Low"),  sum(shp_main$LISA_class_adj=="High-Low")),
  sprintf("| Low-High  | %d | %d |", sum(shp_main$LISA_class_raw=="Low-High"),  sum(shp_main$LISA_class_adj=="Low-High")),
  sprintf("| Non-significant | %d | %d |", sum(shp_main$LISA_class_raw=="Non-significant"), sum(shp_main$LISA_class_adj=="Non-significant")),
  "",
  "---",
  "",
  "## Outputs",
  sprintf("- Figure (PNG+TIFF): `HAV_S2-3t_GiStar_HotspotMap_%s`", TS),
  sprintf("- Figure (PNG+TIFF): `HAV_S2-3hi_LISA_2panel_%s`", TS),
  sprintf("- Source data: `HAV_GiStar_LISA_summary_%s.xlsx`", TS),
  sprintf("- DIR_OUT: `%s`", DIR_OUT),
  "",
  "## Manuscript integration",
  "- Insert as **Supplementary Figure S2-2t** (Gi*) and **Supplementary Figure S2-2hi** (LISA 2-panel)",
  "- Sensitivity 인라인: §2-2.3.5 끝부분에 \"Getis-Ord Gi* and LISA cluster diagnostics confirmed regional hotspot/coldspot patterns consistent with the BYM spatial random effect (Supplementary Figures S2-2t, S2-2hi).\"",
  "",
  "## References (analytical method)",
  "- Getis A, Ord JK. The analysis of spatial association by use of distance statistics. *Geographical Analysis* 1992;24(3):189-206.",
  "- Anselin L. Local indicators of spatial association — LISA. *Geographical Analysis* 1995;27(2):93-115.",
  "",
  "*End of HAV Gi* + LISA Sensitivity Analysis*"
)

writeLines(md, OUT_MD)
cat(sprintf("  ✅ %s\n\n", OUT_MD))

cat(strrep("═", 78), "\n", sep = "")
cat(sprintf("  ✅ HAV Gi* + LISA 분석 완료 (TS: %s)\n", TS))
cat(strrep("═", 78), "\n", sep = "")
}, error=function(e) cat(sprintf("\n\u26a0 [\ud1b5\ud569\ube14\ub85d \uac74\ub108\ub700: GiStar·LISA FigS2] %s\n", conditionMessage(e))))
