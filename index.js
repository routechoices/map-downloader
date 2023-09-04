const app = require("./api");

const port = 3526;

app.listen(port, () => {
    console.log(`App running on http://[::]:${port}`)
})