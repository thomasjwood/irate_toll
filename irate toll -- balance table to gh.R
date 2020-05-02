library(plyr)
library(tidyverse)
library(magrittr)
library(purrrlyr)
library(lsr)
library(broom)
library(emmeans)
library(ggstance)

t1 <- "https://github.com/thomasjwood/irate_toll/raw/master/t1.rds" %>% 
  url %>%
  gzcon %>%
  readRDS

tribble(
  ~cat,   ~value, ~w1, ~w2, ~`0 - 1 corrections`, ~`2 corrections`, ~`3-4 corrections`,
  "total", "",
  t1 %>% 
    use_series(
      "wave"
    ) %>% 
    equals("w1") %>% 
    sum(na.rm = T),
  t1$wave %>% 
    equals("w2") %>% 
    sum,
  t1 %>% 
    filter(
      wave == "w2"
    ) %>% 
    use_series("num_corr") %>% 
    equals(
      t1$num_corr %>% 
        levels %>% 
        extract2(1)
    ) %>% 
    sum(na.rm = T),
  t1 %>% 
    filter(
      wave == "w2"
    ) %>% 
    use_series("num_corr") %>% 
    equals(
      t1$num_corr %>% 
        levels %>% 
        extract2(2)
    )%>% 
    sum(na.rm = T),
  t1 %>% 
    filter(
      wave == "w2"
    ) %>% 
    use_series("num_corr") %>% 
    equals(
      t1$num_corr %>% 
        levels %>% 
        extract2(3)
    ) %>% 
    sum(na.rm = T)
) %>% 
  bind_rows(
    t1 %>% 
      select(
        -num_corr
      ) %>% 
      pivot_longer(
        -c(workerid, wave),
        "cat", "ans"
      ) %>% 
      mutate(
        cat = cat %>% 
          factor(
            cat %>% 
              unique %>% 
              rev
          ),
        value = value %>% 
          factor(
            t1 %>% 
              select(
                age:presvote
              ) %>% 
              map(
                levels
              ) %>% 
              unlist %>% 
              as.character
          )
      ) %>% 
      group_by(
        wave, cat, value
      ) %>% 
      tally %>% 
      spread(
        wave, n
      ) %>% 
      # now we'll left join percentages
      left_join(
        t1 %>% 
          filter(
            wave == "w2"
          ) %>% 
          select(
            -wave
          ) %>% 
          pivot_longer(
            -c(workerid, num_corr),
            "cat", "ans"
          ) %>% 
          na.omit %>%
          mutate(
            cat = cat %>% 
              factor(
                cat %>% 
                  unique %>% 
                  rev
              ),
            value = value %>% 
              factor(
                t1 %>% 
                  select(
                    age:presvote
                  ) %>% 
                  map(
                    levels
                  ) %>% 
                  unlist %>% 
                  as.character
              )
          ) %>% 
          group_by(
            num_corr, cat, value
          ) %>% 
          tally %>% 
          mutate(
            prop = n %>% 
              divide_by(
                n %>% 
                  sum
              ) %>% 
              multiply_by(100) %>% 
              round
          ) %>% 
          select(-n) %>% 
          spread(
            num_corr, prop
          )
      )
  )

# effect sizes

es1 <- t1 %>%
  filter(
    wave == "w2"
  ) %>% 
  select(
    -wave
  ) %>% 
  pivot_longer(
    -c(workerid, num_corr),
    "cat", "ans"
  ) %>% 
  na.omit %>% 
  mutate(
    cat = cat %>% 
      factor(
        cat %>% 
          unique %>% 
          rev
      ),
    value = value %>% 
      factor(
        t1 %>% 
          select(
            age:presvote
          ) %>% 
          map(
            levels
          ) %>% 
          unlist %>% 
          as.character
      )
  ) %>% 
  group_by(
    cat
  ) %>% 
  by_slice(
    function(i){
      
      # i <- es1 %>%
      #   filter(
      #     cat == "presvote"
      #   )
      
      j <- table(factor(i$value), i$num_corr) 
      # 
      # descr::crosstab(
      #   dep = i$value, 
      #   indep = i$num_corr, 
      #   expected = T, 
      #   digits = list(
      #     expected = 0,
      #     others = 3
      #     ),
      #   chisq = T,
      #   prop.chisq = T
      #   )
      
      str_c(
        j %>% 
          chisq.test %>% 
          tidy %>%
          use_series(statistic) %>% 
          round(2),
        j %>% 
          chisq.test %>% 
          tidy %>% 
          use_series(p.value) %>% 
          gtools::stars.pval() %>% 
          str_trim,
        " (",
        j %>% 
          cramersV %>% 
          round(2) %>% 
          as.character %>% 
          str_replace(fixed("0."), "."),
        ")"
      )
    }, 
    .collate = "rows"
  )

# figure 2


t2 <- "https://github.com/thomasjwood/irate_toll/raw/master/t2.rds" %>% 
  url %>%
  gzcon %>% 
  readRDS

l2 <- t2 %>% 
  group_by(issue) %>% 
  by_slice(
    
    function(i)
      str_c(
        "ans_num ~ cond",
        c("",
          "*ideo3")
      ) %>%
      map(
        function(j)
          lm(j, i)
      ) %>% 
      map2(
        c(pairwise ~ cond,
          pairwise ~ cond | ideo3),
        function(k, l)
          emmeans(
            k, l, data = i
          )
      )
  )

bt1 <- l2$issue %>% 
  mapvalues(
    l2$issue %>% 
      unique,
    c("China holds the majority of U.S. debt",
      "The trade deficit with China is a sum of money we pay to China, for which we receive nothing in return",
      "The U.S. manufacturing sector has been steadily losing jobs over the last decade",
      "US workers' pay is at an all-time low")
  ) %>% 
  map2_dfr(
    l2$.out, 
    function(h, i)
      
      i %>%
      map_dfr(
        ~extract2(., 1) %>% 
          tidy
      ) %>% 
      mutate(
        ideo3 = ideo3 %>% 
          is.na %>% 
          ifelse(
            "Overall",
            ideo3 %>% 
              as.character
          ) %>% 
          factor(
            c("Overall",
              "Conservative",
              "Moderate",
              "Liberal")
          ),
        issue = h
      )
  ) %>% 
  mutate(
    issue = issue %<>%
      fct_reorder(
        estimate
      )
  )

bt1 %>% 
  ggplot() +
  geom_linerangeh(
    aes(
      xmin = conf.low, xmax = conf.high, 
      y = ideo3, group = cond,
      linetype = cond
    ),
    position = position_dodgev(.2)
  ) +
  geom_point(
    aes(
      x = estimate,
      y = ideo3,
      group = cond
    ),
    shape = 21,
    fill = "white",
    size = 6,
    position = position_dodgev(.2)
  ) +
  geom_point(
    aes(
      x = estimate,
      y = ideo3,
      group = cond
    ),
    shape = 21,
    fill = "white",
    color = "white",
    size = 5,
    position = position_dodgev(.2)
  ) +
  geom_text(
    aes(
      x = estimate,
      y = ideo3,
      group = cond,
      label = estimate %>% 
        round(1)
    ),
    size = 2,
    position = position_dodgev(.2),
    family = "Roboto"
  ) +
  geom_text(
    aes(
      x, y, label = label
    ), 
    size = 3,
    nudge_y = .65,
    data = data.frame(
      x = c(2.1, 3.4),
      y = c("Liberal", "Liberal"),
      label = c(
        "Correction",
        "No Correction"
      )
    ) %>% 
      mutate(
        issue = bt1$issue %>% 
          levels %>% 
          extract2(1) %>% 
          factor(
            bt1$issue %>% 
              levels
          )
      ),
    fontface = "bold.italic",
    family = "Roboto",
  ) +
  facet_grid(
    issue ~ ., 
    labeller = label_wrap_gen(width = 30),
    space = "free_y"
  ) +
  labs(
    x = "Agreement (7pt scale)",
    y = ""
  ) 

# correction effects

l4_1 <- "https://github.com/thomasjwood/irate_toll/raw/master/t3_1.rds" %>% 
  url %>% 
  gzcon %>% 
  readRDS %>% 
  by_slice(
    function(i)
      
      c("ans ~ num_corr",
        "ans ~ ideo3 * num_corr") %>% 
      map(
        function(j)
          lm(j, i)
      ) %>% 
      map2(
        c(revpairwise ~ num_corr,
          revpairwise ~ num_corr | ideo3),
        function(k, l)
          emmeans(
            k, l, 
            data = i,
            at = list(
              num_corr = c(0, 4)
            )
          )
      )
  ) 

bt4 <- l4_1$dv %>% 
  map2_dfr(
    l4_1$.out,
    function(i, j)
      
      # i <- l4$dv[[1]]
      # j <- l4$.out[[1]]
      
      j %>% 
      map_dfr(
        function(k)
          k %>% 
          extract2(1) %>% 
          tidy
      ) %>% 
      mutate(
        ideo3 = ideo3 %>% 
          is.na %>% 
          ifelse(
            "Overall",
            ideo3 %>% 
              as.character
          ) %>% 
          factor(
            c("Overall",
              "Conservative",
              "Moderate",
              "Liberal")
          ),
        dv = i
      ) %>%
      left_join(
        j %>% 
          map_dfr(
            function(k)
              k %>% 
              extract2(2) %>% 
              tidy
          ) %>% 
          select(
            level2:std.error, p.value, ideo3
          ) %>% 
          mutate(
            ideo3 = ideo3 %>% 
              is.na %>% 
              ifelse(
                "Overall",
                ideo3 %>% 
                  as.character
              ) %>% 
              factor(
                c("Overall",
                  "Conservative",
                  "Moderate",
                  "Liberal")
              ),
            level2 = level2 %>% 
              as.numeric
          ) %>% 
          rename(
            num_corr = level2, 
            diff_est = estimate,
            diff_se = std.error
          )
      )
  ) %>% 
  mutate(
    dv = dv %>% 
      mapvalues(
        c("overall", "pol_agree_good", "pol_create_jobs", "pol_economy_grow", 
          "pol_family_help", "pol_goods_price", "pol_worker_wages"),
        c("overall",
          "...between US and other countries a good thing",
          "...create jobs",
          "...make American economy grow",
          "...help my family's financial situation",
          "Free trade agreements...reduce prices for products sold in US",
          "...increase wages for American workers")
      ),
    dv = dv %>% 
      fct_reorder(
        estimate, 
        .desc = T
      ),
    num_corr = num_corr %>% 
      factor
  )

bt4 %<>% 
  left_join(
    bt4 %>% 
      select(
        estimate, dv, ideo3
      ) %>% 
      group_by(
        dv, ideo3
      ) %>% 
      summarize(
        corr_mid = estimate %>% 
          mean,
        num_corr = c("0") %>% 
          factor(
            c("0", "4")
          )
      )
  )


bt4$p.value[seq(2, 56, 2)] = bt4$p.value[seq(1, 55, 2)]

bt4$dv %<>% 
  mapvalues(
    "overall",
    "Overall"
  ) %>% 
  factor(
    bt4$dv %>% 
      levels %>% 
      extract(-4) %>% 
      c("Overall")
  )

bt4_1 <- bt4 %>%
  filter(
    p.value <= .05
  ) %>% 
  group_by(
    ideo3,
    dv
  ) %>% 
  summarize(
    start = estimate %>% 
      min,
    end = estimate %>% 
      max
  )

bt4 %>% 
  ggplot() +
  geom_linerangeh(
    aes(
      xmin =  start,
      xmax = end,
      y = ideo3 %>% 
        as.numeric %>% 
        add(.2)
    ),
    data = bt4_1,
    color = "grey50", 
    size = .4,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(x = start,
        xend = start,
        y = ideo3 %>%
          as.numeric %>% 
          subtract(.05),
        yend = ideo3 %>%
          as.numeric %>% 
          add(.2)
    ),
    data = bt4_1,
    size = .2,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(x = end,
        xend = end,
        y = ideo3 %>%
          as.numeric %>% 
          add(.05),
        yend = ideo3 %>%
          as.numeric %>% 
          add(.2)
    ),
    data = bt4_1,
    size = .2,
    linetype = "dashed"
  ) +
  geom_linerangeh(
    aes(xmin = conf.low, 
        xmax = conf.high,
        y = ideo3 %>% 
          as.numeric,
        linetype = num_corr),
    position = position_dodgev(.2)
  ) +
  geom_label(
    aes(
      corr_mid, ideo3 %>% 
        as.numeric,
      label = diff_est %>% 
        round(2) %>% 
        as.character %>% 
        str_sub(2) %>% 
        str_pad(
          width = 3,
          side = "right",
          pad = "0"
        ) %>% 
        str_c(
          c("***", 
            "**",
            "*",
            "") %>% 
            extract(
              p.value %>% 
                findInterval(
                  c(-Inf, .001, .01, .05, Inf)
                )
            )
        )
    ),
    fill = "grey95",
    label.size = 0, 
    nudge_y = .2,
    size = 2.25,
    family = "Roboto",
    data = bt4 %>% 
      filter(
        p.value <= .05
      )
  ) +
  geom_point(
    aes(
      estimate, ideo3%>% 
        as.numeric,
      group = num_corr
    ),
    position = position_dodgev(.2),
    shape = 21,
    fill = "white",
    size = 4.5
  ) +
  geom_point(
    aes(
      estimate, ideo3 %>% 
        as.numeric,
      group = num_corr
    ),
    position = position_dodgev(.2),
    shape = 21,
    fill = "white",
    color = "white",
    size = 3.75
  ) +
  geom_text(
    aes(
      estimate, 
      ideo3 %>% 
        as.numeric, 
      label = estimate %>% 
        round(2) %>% 
        as.character %>% 
        str_sub(2) %>% 
        str_pad(
          width = 3,
          side = "right",
          pad = "0"
        ),
      group = num_corr
    ),
    family = "Roboto",
    position = position_dodgev(.2),
    size = 1.5
  ) +
  geom_text(
    aes(
      x, y, label = label
    ),
    data = data.frame(
      x = c(.48, .78),
      y = c(2, 2) %>% 
        add(
          c(-.25, -.15)
        ),
      label = c("0\ncorrections",
                "4\ncorrections")
    ) %>% 
      mutate(
        dv = bt4$dv %>% 
          levels %>% 
          extract2(1) %>% 
          factor(
            bt4$dv %>% 
              levels 
          )
      ), 
    lineheight = .7,
    size = 2.5,
    fontface = "bold.italic",
    family = "Roboto",
  ) +
  facet_grid(
    . ~ dv,
    space = "free_x",
    scales = "free_x",
    labeller = label_wrap_gen(width = 15)
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, .2),
    labels = seq(0, 1, .2) %>% 
      as.character %>% 
      str_replace_all(fixed("0."), ".")
  ) +
  scale_y_continuous(
    breaks = 1:4,
    labels = bt4$ideo3 %>% 
      levels
  ) +
  labs(
    x = "Probability of agreeing with statement",
    y = ""
  ) +
  theme_minimal(base_family = "Roboto") +
  theme(
    panel.background  = 
      element_rect(color = "grey95",
                   fill = "grey95"),
    panel.grid = element_blank(),
    strip.text.x = element_text(
      angle = 0,
      face = bt4$dv %>% 
        levels %>% 
        equals("Overall") %>% 
        ifelse(
          "bold", "plain"
        )
    ),
    strip.background = element_rect(color = "grey95",
                                    fill = "grey95"),
    legend.background = element_rect(color = "white",
                                     fill = "white"),
    legend.margin = margin(0, 0, 0, 0, "cm"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.caption = element_text(face = "italic")
  )
