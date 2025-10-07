library(shiny)
runApp()

library(rsconnect)
rsconnect::setAccountInfo(name='mouselivermultiomics',
                          token='3518F9D3CF80607D55A2AB67DDDF3B3A',
                          secret='QCDRkc4QXhqsTHYACewpyeezcqI+C8ZCN6+Ihyi+')
deployApp()
