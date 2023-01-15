library(usethis)

use_git_config(user.name="Junier Marrero Gutiérrez",
               user.email="juniermarreroguti@gmail.com")


usethis::create_github_token()

# Token
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

#
gitcreds::gitcreds_set()


# Reinicie o R usando o RStudio: CTRL + SHIFT + F10
usethis::git_sitrep()

